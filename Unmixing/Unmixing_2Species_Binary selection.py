
"""
 Script to unmix two species in confocal-FLIM images using phasor analysis.
 This script allows the user to select folders containing control images for two fluorophores and a folder with mixed images.
 It calculates the phasor coordinates, performs unmixing, and generates output images and plots.

  The script outputs:
    - the phasor distribution of the mixed images colored by the unmixing results
    - the mixed images colored by the unmixing results and the separated fractions
    - the phasor distribution of the control images
    - Composite image of the unmixed fluorophores with the third fraction removed
"""

import os
import glob
import numpy
import easygui
import pandas as pd
import scipy
import time
import matplotlib
import matplotlib.pyplot as plt
import skimage
import tifffile
from sklearn.cluster import KMeans
from matplotlib.patches import Circle, Ellipse
#plt.style.use('dark_background')
from sys import path as path1;
Functionspath=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Functions")
path1.append(Functionspath)
from Main_functions import (line_equation, to_polar_coord, polar_to_cart, load_image,select_channel, get_foreground)
from Phasor_functions import Median_Phasor,unmix2species
from tiffwrapper import imsave,LifetimeOverlayer


#plt.style.use('dark_background')
# ------------------ Default Input variables ----------------
params_dict = {
    "smooth_factor": 0.2,  # % 'Smoothing factor:'
    "phasor_smooth_cycles": 1,  # % 'Smoothing cycles phasor:'
    "foreground_threshold": 10,
    "harm1": 1,
}

matplotlib.rcParams['axes.linewidth'] = 0.8
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    name="coolspring",
    colors=["#00ffffff", "#ff00ffff"]
)
matplotlib.colormaps.register(cmap=cmap, force=True)
matplotlib.colormaps.register(cmap=cmap.reversed(), force=True)
matplotlib.rcParams['axes.linewidth'] = 0.8

# -----------------------------------------------------------

# Select the folders containing the control images and the mixed images. Open a dialog box to select the folders.
#f1=easygui.diropenbox(default=os.path.expanduser("~Desktop"), title="Select folder containing the control images for the first fluorophore")
#f2=easygui.diropenbox(default=os.path.expanduser("~Desktop"), title="Select folder containing the control images for the second fluorophore")
#f3=easygui.diropenbox(default=os.path.expanduser("~Desktop"), title="Select folder containing the mixed images (the mixture of the two fluorophores)")
f1="D:\\Zenodo\\FarRedDyes_SingleLabel\\FixedCells\\alphaTubulin_Alexa647\\"
f2="D:\\Zenodo\\FarRedDyes_SingleLabel\\FixedCells\\rabBassoon STAR635P\\"
f3="D:\\Zenodo\\FarRedDyes_DualLabel\\FixedCells\\Tubulin_Alexa647_Bassoon_STAR635P\\"

filenamescontrol = [f1, f2]
filenamemixed=f3

# Labels for the legend
# These labels correspond to the fluorophores used in the control images and the mixed images.
colors=['magenta',  'c' ,'lightgreen']
cmapbin = matplotlib.colors.ListedColormap(["magenta","cyan","gray"])
normmpl = matplotlib.colors.BoundaryNorm([0,1,2], cmapbin.N)
labels = ['Bassoon CF594', 'PSD95 STAR Orange', 'Mixture']

# Which channels to use (keys in the msr files, channel numbers in the tiff files)
keys = ['Confocal_561 {11}', 'Confocal_561 {11}','Confocal_561 {11}']

#keys = [ 'Conf_635P {2}','Conf_635P {2}','Conf_635P {2}' ]
keys = [0,0,0]# For Tiff files
msrfiles = []

# -----------------------------------------------------------
# Ask user for name of output folder and create it on the Desktop
savefolder=str(input("Name of Output folder: "))
savefolder = os.path.join(os.path.expanduser("~/Desktop"), "Unmixing_"+savefolder+"_Confocal")
os.makedirs(savefolder, exist_ok=True)

# Make list of all the images in the control folders and ask user to select the image to use as control
for filename in filenamescontrol :
   
    extension = ".msr"
    path=os.path.join(filename,"*.msr")
    images = glob.glob(path)
    print('There are ',len(images), ' msr files in this folder')
    if len(images) == 0:
        path=os.path.join(filename,"*.tiff")
        images = glob.glob(path)
        print('There are ',len(images), ' tiff files in this folder')
        extension = ".tiff"
    for idx, imagei in enumerate(images):
        print(idx,os.path.basename(imagei)) 
    numim = int(input('Please enter the image number (1st=0): '))
    image = images[numim]
    msrfiles.append(image)
print(msrfiles)

path = os.path.join(filenamemixed, '*'+extension)
mixedimages = glob.glob(path)

# Create Figure and axes for the phasor plot
fig4,ax_scatter = plt.subplots(figsize=(2,2))
ax_scatter.set_xlim(0, 1.05)
ax_scatter.set_ylim(-0.05, 1)

ax_scatter.set_xlabel('g')
ax_scatter.set_ylabel('s')
# Plot the universal semi-circle in the phasor plot
edge = numpy.linspace(start=0, stop=15, num=200)
theta = numpy.linspace(0, numpy.pi, 100)
r = 0.5
x1 = r*numpy.cos(theta) + 0.5
x2 = r*numpy.sin(theta)
ax_scatter.plot(x1,x2, color = "black", ls = "--",linewidth=0.8)
CoM_x, CoM_y = [], []
d_melange = pd.DataFrame(columns=['g', 's'])
 # Create a legend file 
with open(os.path.join(savefolder,'legend.txt'),'w') as data:
    data.write("Controls\n")

scatterlist = []
for i, msr in enumerate(msrfiles) : 
    df = pd.DataFrame(columns=['x','y'])
    dg = pd.DataFrame(columns=['g', 's'])
     # Write info about the control image used to to the legend file
    with open(os.path.join(savefolder,'legend.txt'),'a') as data:
        data.write("{}\t{}\t{}\n".format(labels[i],keys[i],msr))

    # Load the control image and select the channel  
    imagemsr=load_image(msr)
    image1=select_channel(imagemsr,keys[i])
  
    print(image1.shape)
    imsum = image1.sum(axis=2)
    imsum = imsum.astype('int16')
    
    seuil = get_foreground(imsum)
    print("Caclulation for an image of shape", image1.shape, "...")
    params_dict["foreground_threshold"] = seuil
    print("foreground_threshold=", params_dict["foreground_threshold"])
    # Calulate the phasor coordinates for the foreground of the control image
    x,y,g_smoothed,s_smoothed, original_idxes= Median_Phasor(image1, params_dict, **params_dict)
    df['x']=x.flatten()
    df['y']=y.flatten()
    m, phi = to_polar_coord(df['x'], df['y'])
    g,s =polar_to_cart(m, phi)
    dg['g'], dg['s'] = g, s

    # Calculate the centroid of the phasor distribution using KMeans clustering
    kmeans = KMeans(n_clusters = 1, init = 'k-means++', random_state = 42)
    y_kmeans = kmeans.fit_predict(dg)
    CoM_x.extend(kmeans.cluster_centers_[:, 0][:].tolist())
    CoM_y.extend(kmeans.cluster_centers_[:, 1][:].tolist())

    # Plot the phasor distribution for the control image
    a=ax_scatter.scatter(g, s, s=1, c=colors[i], alpha=0.1,label=labels[i],rasterized=True)
    scatterlist.append(a)

xaxis = numpy.linspace(0, 1, 100)

## Calculate the nearest point on the semi-circle to the centroids of the phasor distributions to define the points Pn and P2
norm = numpy.sqrt((CoM_x[0] - 0.5) ** 2 + (CoM_y[0] ** 2))
Pn_x = 0.5 + (r * (CoM_x[0] - 0.5) / norm)
Pn_y = 0 + r * (CoM_y[0] - 0) / norm
P_n = numpy.array([Pn_x, Pn_y])
norm = numpy.sqrt((CoM_x[1] - 0.5) ** 2 + (CoM_y[1] ** 2))
P2_x = 0.5 + (r * (CoM_x[1] - 0.5) / norm)
P2_y = 0 + r * (CoM_y[1] - 0) / norm
p2 = numpy.array([P2_x, P2_y])

dim_cercle = {
    "x0" : P2_x, "y0" : P2_y, "r"  : 0.2,
    "x1" : Pn_x, "y1" : Pn_y, "r1" : 0.125,

}
pnscatter=ax_scatter.scatter(Pn_x,Pn_y, c='darkred' )
p2scatter=ax_scatter.scatter(P2_x,P2_y, c='darkred' )
circle1 = Circle((dim_cercle["x0"], dim_cercle["y0"]), dim_cercle["r"], fill=False, color= 'yellow')
ax_scatter.add_patch(circle1)
circle2 = Circle((dim_cercle["x1"], dim_cercle["y1"]), dim_cercle["r1"], fill=False, color = 'blue')
ax_scatter.add_patch(circle2)
fig4.savefig(os.path.join(savefolder,
                            "Phasor_2species_Confocal_ControlsOnly.pdf"),
                transparent='True',
                bbox_inches="tight", dpi=900)
 # Clear the scatter points and lines from the phasor plot 
pnscatter.remove()
p2scatter.remove()
ax_scatter.lines[-1].remove()
t = [scatter.remove() for scatter in scatterlist]

 # Loop over the mixed images and perform unmixing
for m,mixedimage in enumerate(mixedimages):
    print("***********************************************************")
    print("Working on image number ",m," out of ",len(mixedimages))
    print("***********************************************************")
    d_melange = pd.DataFrame(columns=['g', 's'])
    df = pd.DataFrame(columns=['x', 'y'])
    dg = pd.DataFrame(columns=['g', 's'])

    # Load the mixed image and select the channel
    imagemsr = load_image(mixedimage)
    image1 = select_channel(imagemsr, keys[2])
    print(image1.shape)

    imsum = image1.sum(axis=2)
    imsum = imsum.astype('int16')

    params_dict["foreground_threshold"] = 5
    
   
    print("foreground_threshold=", params_dict["foreground_threshold"])
    # Calculate the phasor distribution of the foreground for the mixed image
    x, y, g_smoothed, s_smoothed, original_idxes = Median_Phasor(image1, params_dict, **params_dict)
    df['x'] = x.flatten()
    df['y'] = y.flatten()

    # Calibrate the phasor distribution using the IRF
    m, phi = to_polar_coord(df['x'], df['y'])
    g, s = polar_to_cart(m, phi)
    dg['g'], dg['s'] = g, s
    d_melange['g'], d_melange['s'] = g, s

    # Calculate the centroid of the phasor distribution using KMeans clustering
    kmeans = KMeans(n_clusters=1, init='k-means++', random_state=42)
    y_kmeans = kmeans.fit_predict(dg)
    CoM_x.extend(kmeans.cluster_centers_[:, 0][:].tolist())
    CoM_y.extend(kmeans.cluster_centers_[:, 1][:].tolist())

    p3 = d_melange[['g', 's']].to_numpy()
    Solve=numpy.zeros(p3.shape[0])
    print("P3 shape:", p3.shape)

    xyTrue1 = ((p3[:, 0] - dim_cercle["x0"]) ** 2 + (p3[:, 1] - dim_cercle["y0"]) ** 2) < dim_cercle["r"]**2
    xyTrue2 = ((p3[:, 0] - dim_cercle["x1"]) ** 2 + (p3[:, 1] - dim_cercle["y1"]) ** 2) < dim_cercle["r1"]**2
    Solve[xyTrue1==True]=0
    Solve[xyTrue2==True]=1
    Solve[(xyTrue1==False) & (xyTrue2==False)]=2
    print("Number of pixels assigned to species 1:", numpy.sum(Solve==0))
    print("Number of pixels assigned to species 2:", numpy.sum(Solve==1))
    print("Number of pixels assigned to neither:", numpy.sum(Solve==2))




    # Plot the phasor distribution of the mixed image colored by the unmixing results and save the figure
    mixphasor = ax_scatter.scatter(g, s, s=1, c=Solve,cmap=cmapbin, alpha=0.1,label="Mixture",rasterized=True)
    
    pnscatter = ax_scatter.scatter(Pn_x, Pn_y, c='darkred')
    p2scatter = ax_scatter.scatter(P2_x, P2_y, c='darkred')
    ax_scatter.plot(x1,x2, color = "black", ls = "--",linewidth=0.8)

    fig4.savefig(os.path.join(savefolder, "Phasor_2species_Confocal_{}.pdf".format(os.path.basename(mixedimage).split(".msr")[0])), transparent='True',
                    bbox_inches="tight",dpi=900)
    
    # Clear the scatter points and lines from the phasor plot
    pnscatter.remove()
    p2scatter.remove()

    mixphasor.remove()
    fraction1 = numpy.array(imsum.flatten()*0)
    fraction2 =numpy.array(imsum.flatten()*0)
    print("Original idxes length:", len(original_idxes))
    for idx,j in enumerate(original_idxes) :
        if Solve[idx]==0:
            fraction1[j] = 1
        elif Solve[idx]==1:
            fraction2[j] = 1
    fraction1 = fraction1.reshape(imsum.shape)
    fraction2 = fraction2.reshape(imsum.shape)



    # Create a figure to display the images and the unmixing results        
    fig_im3, ax_im3 = plt.subplots(ncols=4,nrows=1,figsize=(12,8))
    ax_im3[0].axis('off')
    ax_im3[1].axis('off')
    ax_im3[2].axis('off')
    ax_im3[3].axis('off')

    imsum_flat3 =ax_im3[1].imshow(fraction1*imsum, cmap='hot')
    cbar3 =fig_im3.colorbar(imsum_flat3, ax=ax_im3[1],fraction=0.05, pad=0.01)
    imsum_flat5 =ax_im3[2].imshow(fraction2*imsum, cmap='hot')
    cbar2 =fig_im3.colorbar(imsum_flat5, ax=ax_im3[2],fraction=0.05, pad=0.01)
    imsum_flat6 =ax_im3[0].imshow(imsum, cmap='hot')
    cbar1 =fig_im3.colorbar(imsum_flat6, ax=ax_im3[0],fraction=0.05, pad=0.01)
    overlayer = LifetimeOverlayer(fraction1, imsum/imsum.max(), cname='coolspring')
    lifetime_rgb, cmap = overlayer.get_overlay(
        lifetime_minmax=(0., 1),
        intensity_minmax=(0, 0.5) # inTensity saturated to get more bright regions
                )
    imsum_flat5 =ax_im3[3].imshow(lifetime_rgb)
    cbar =fig_im3.colorbar(cmap, ax=ax_im3[3],fraction=0.05, pad=0.01)


    imagecomp=numpy.dstack((fraction1*imsum, fraction2*imsum))
    imagecomp=numpy.moveaxis(imagecomp,2,0)
    filenameout =  os.path.join(savefolder,os.path.basename(mixedimage).split(extension)[0] + "_ConfCentroidsCircle_UnmixedComposite.tiff")
    imsave(file=filenameout, data=imagecomp.astype(numpy.uint16), composite=True, luts=("Cyan Hot","Magenta Hot"), pixelsize=(20E-3,20E-3))

    filenameout = os.path.join(savefolder,os.path.basename(mixedimage).split(extension)[0] + "_ConfCentroidsCircle_MixedIntensity.tiff")
    print(filenameout)
    imsave(file=filenameout, data=imsum.astype(numpy.uint16), luts="Red Hot", pixelsize=(20E-3,20E-3))


    filenameout =  os.path.join(savefolder,os.path.basename(mixedimage).split(extension)[0] + "_ConfCentroidsCircle_f1f2.tiff")
    imagecomp=numpy.dstack((fraction2,fraction1))
    imagecomp=numpy.moveaxis(imagecomp,2,0)
    tifffile.imwrite(filenameout, imagecomp)

    filenameout = os.path.join(savefolder,
                                os.path.basename(mixedimage).split(extension)[0] + "_Conf2species_LineControls_F1Overlay.tiff")

    tifffile.imwrite(filenameout, lifetime_rgb.astype(numpy.float32))




    fig_im3.savefig(os.path.join(savefolder,'Images__SeparateConf_CentroidsCircle'+os.path.basename(mixedimage).split(extension)[0] +'.pdf'),transparent='True', bbox_inches="tight")


