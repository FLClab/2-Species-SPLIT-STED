
"""
 Script to unmix two species in STED-FLIM images using phasor analysis.
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
import time
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import easygui
import seaborn
import tifffile
from sklearn.cluster import KMeans
import matplotlib
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
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    name="coolspring",
    colors=["#00ffffff", "#ff00ffff"]
)
matplotlib.colormaps.register(cmap=cmap, force=True)
matplotlib.colormaps.register(cmap=cmap.reversed(), force=True)
matplotlib.rcParams['axes.linewidth'] = 0.8

# -----------------------------------------------------------

# Select the folders containing the control images and the mixed images. Open a dialog box to select the folders.
f1=easygui.diropenbox(default=os.path.expanduser("~Desktop"), title="Select folder containing the control images for the first fluorophore")
f2=easygui.diropenbox(default=os.path.expanduser("~Desktop"), title="Select folder containing the control images for the second fluorophore")
f3=easygui.diropenbox(default=os.path.expanduser("~Desktop"), title="Select folder containing the mixed images (the mixture of the two fluorophores)")
filenamescontrol = [f1, f2]
filenamemixed=f3



# Labels for the legend
# These labels correspond to the fluorophores used in the control images and the mixed images.
labels = ['Bassoon CF594', 'PSD95 STOrange', 'Mixture']
colors=['magenta',  'c' ,'springgreen']

# Which channels to use (keys in the msr files, channel numbers in the tiff files)

#keys=[1,1,1] # For Tiff files

#keys=['STED640 {10}', 'STED640 {10}', 'STED640 {10}']
#keys = ['STED_635P {2}', 'STED_635P {2}', 'STED_635P {2}']
keys = ['STED 561 {11}', 'STED 561 {11}','STED 561 {11}']


# List of STED depletion powers to analyze and the ID of the images to use as controls.
powers=[[10,[0,0]],[20,[0,0]],[30,[0,0]],[40,[0,0]]]
#powers=[[5,[0,0]],[10,[0,0]],[15,[0,0]],[20,[0,0]]]

# Aske the user to input the output folder name. Create the output folder on the Desktop.
savefolder=str(input("Name of Output folder: "))
savefoldermain = os.path.join(os.path.expanduser("~/Desktop"), "Unmixing_"+savefolder+"_2Species")
os.makedirs(savefoldermain, exist_ok=True)

#Make figure with colorbar for phasor legend
img = plt.imshow(numpy.array([[0,1]]), cmap="cool")
img.set_visible(False)

plt.colorbar(orientation="vertical")
plt.savefig(os.path.join(savefoldermain,"Legend.pdf"),transparent=True)

# Loop over the powers and process the images for each power.
for power in powers:
    print(power)

    STEDPOWER=power[0]
    # Create the output folder for the current power.
    savefolder = os.path.join( savefoldermain,"{}PercentSTED".format(STEDPOWER))
    os.makedirs(savefolder, exist_ok=True)


    msrfiles = []
    for k,filename in enumerate(filenamescontrol) :
        # Make list of all the images in the folder
        extension = ".msr"
        path = os.path.join(filename, '*_{}PercentSTED.msr'.format(STEDPOWER) )
        images = glob.glob(path)
        print('There are ',len(images), ' msr files in this folder')
        if len(images) == 0:
            path = os.path.join(filename, '*_{}PercentSTED.tiff'.format(STEDPOWER) )
            images = glob.glob(path)
            print('There are ',len(images), ' tiff files in this folder')
            extension = ".tiff"
        for i,imagei in enumerate(images):
            print(i,os.path.basename(imagei))
        numim = power[1][k]
        image = images[numim]
        msrfiles.append(image)
    print(msrfiles)

    path = os.path.join(filenamemixed, '*_{}PercentSTED{}'.format(STEDPOWER,extension))
    mixedimages = glob.glob(path)

# Create Figure and axes for the phasor plot
    fig, ax = plt.subplots()
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
    ax_scatter.plot(x1,x2, color = "k", ls = "--",linewidth=0.8)
    

    # Create a legend file 
    with open(os.path.join(savefolder,'legend.txt'),'w') as data:
        data.write("Controls\n")

    scatterlist = []
    CoM_x, CoM_y = [], []
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
        a=ax_scatter.scatter(g, s, s=1, c=colors[i], alpha=0.1, label=labels[i],rasterized=True)
        scatterlist.append(a)


    ## Use the centroids of the phasor distributions to define the points Pn and P2
    Pn_x, Pn_y = CoM_x[0], CoM_y[0]
    P_n =numpy.array([Pn_x, Pn_y])

    P2_x =CoM_x[1]
    P2_y =CoM_y[1]
    p2=numpy.array([P2_x,P2_y])

    ## Create linear trajectory connecting Pn and p2
    xaxis = numpy.linspace(0, 1, 100)
    x2, y2 = [Pn_x, P2_x], [Pn_y,P2_y]
    m2, c2 = line_equation(x2, y2)
    y2 = m2 * xaxis + c2

    # Plot the linear trajectory in the phasor plot and save the figure
    p2pnline = ax_scatter.plot([Pn_x, P2_x], [Pn_y, P2_y], c='dodgerblue')
    pnscatter = ax_scatter.scatter(Pn_x, Pn_y, s=50, c='gold')
    p2scatter = ax_scatter.scatter(P2_x, P2_y, s=50, c='gold')
    fig4.savefig(os.path.join(savefolder, "Phasor_2species_ControlsOnly.pdf"), transparent='True',
                    bbox_inches="tight",dpi=900)
    
    # Clear the scatter points and lines from the phasor plot 
    t = [scatter.remove() for scatter in scatterlist]
    lines = [mpatches.Patch(color=colors[j], label=labels[j]) for j in range(len(labels))]
    pnscatter.remove()
    p2scatter.remove()
    ax_scatter.lines[-1].remove()


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
        imsum = image1[:,:,10:111].sum(axis=2)
        imsum = imsum.astype('int16')

        # Calculate the phasor distribution of the foreground for the mixed image
        params_dict["foreground_threshold"] = 5

        print("foreground_threshold=", params_dict["foreground_threshold"])

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
        # Perform unmixing using the phasor coordinates and the points Pn and P2
        imsum_flat_lin1,imsum_flat_lin2,Solve=unmix2species(p3, original_idxes,image1,P_n,p2)
        fraction1 = imsum_flat_lin2.copy()
        fraction2 = imsum_flat_lin1.copy()
        imsum_flat_lin1*= imsum
        imsum_flat_lin2*= imsum

        # Plot the phasor distribution of the mixed image colored by the unmixing results and save the figure
        mixphasor = ax_scatter.scatter(g, s, s=1,c=Solve[1,:],cmap="cool",rasterized=True,label="Mixture")

        lines = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(len(labels))]
        p2pnline = ax_scatter.plot([Pn_x, P2_x], [Pn_y, P2_y], c='dodgerblue')
        pnscatter = ax_scatter.scatter(Pn_x, Pn_y, s=50, c='gold')
        p2scatter = ax_scatter.scatter(P2_x, P2_y, s=50, c='gold')
        fig4.savefig(os.path.join(savefolder, "Phasor_2species_{}.pdf".format(os.path.basename(mixedimage).split(extension)[0])), transparent='True',
                        bbox_inches="tight",dpi=900)
        
        # Clear the scatter points and lines from the phasor plot
        pnscatter.remove()
        p2scatter.remove()
        ax_scatter.lines[-1].remove()
        mixphasor.remove()

        # Create a figure to display the images and the unmixing results
        fig_im3, ax_im3 = plt.subplots(ncols=4,nrows=1,figsize=(12,8))
        ax_im3[0].axis('off')
        ax_im3[1].axis('off')
        ax_im3[2].axis('off')
        ax_im3[3].axis('off')

        imsum_flat3 =ax_im3[1].imshow(imsum_flat_lin1, cmap='hot')
        cbar3 =fig_im3.colorbar(imsum_flat3, ax=ax_im3[1],fraction=0.05, pad=0.01)
        imsum_flat5 =ax_im3[2].imshow(imsum_flat_lin2, cmap='hot')
        cbar2 =fig_im3.colorbar(imsum_flat5, ax=ax_im3[2],fraction=0.05, pad=0.01)
        imsum_flat6 =ax_im3[0].imshow(imsum, cmap='hot')
        cbar1 =fig_im3.colorbar(imsum_flat6, ax=ax_im3[0],fraction=0.05, pad=0.01)

        
        overlayer = LifetimeOverlayer(fraction1, imsum/imsum.max(),  cname="coolspring")
        lifetime_rgb, cmap = overlayer.get_overlay(
            lifetime_minmax=(0., 1),
            intensity_minmax=(0, 0.5) # inTensity saturated to get more bright regions
                    )
        imsum_flat5 =ax_im3[3].imshow(lifetime_rgb)
        cbar =fig_im3.colorbar(cmap, ax=ax_im3[3],fraction=0.05, pad=0.01)
        filenameout = os.path.join(savefolder,
                                    os.path.basename(mixedimage).split(extension)[0] + "_STED2species_F1Overlay.tiff")

        tifffile.imwrite(filenameout, lifetime_rgb.astype(numpy.float32))
    
        lifetime=numpy.dstack(( fraction2, fraction1, fraction1*0))


        red=lifetime[:,:,1]+lifetime[:,:,2]
        green=lifetime[:,:,0]+lifetime[:,:,2]
        blue=lifetime[:,:,0]+lifetime[:,:,1]
        
        lifetime_rgb =numpy.dstack((red,green,blue))
        #lifetime_rgb=numpy.where(lifetime_rgb>1, 1,lifetime_rgb)
        print(lifetime_rgb.shape,lifetime_rgb.dtype)
        print(numpy.min(lifetime_rgb[:,:,0]),numpy.max(lifetime_rgb[:,:,0]))  
        print(numpy.min(lifetime_rgb[:,:,1]),numpy.max(lifetime_rgb[:,:,1]))  
        print(numpy.min(lifetime_rgb[:,:,2]),numpy.max(lifetime_rgb[:,:,2]))  

        plt.figure()
        plt.imshow(lifetime_rgb)
        plt.axis('off')

        plt.savefig(os.path.join(savefolder, os.path.basename(mixedimage).split(extension)[0] + "_STED2species_lifetimergb.pdf"), transparent='True', bbox_inches="tight")
        


        imagecomp=numpy.dstack((imsum_flat_lin2,imsum_flat_lin1))

        imagecomp=numpy.moveaxis(imagecomp,2,0)

        filenameout =  os.path.join(savefolder,os.path.basename(mixedimage).split(extension)[0] + "_STED2species_UnmixedComposite.tiff")
        imsave(file=filenameout, data=imagecomp.astype(numpy.uint16), composite=True, luts=("magenta","cyan"), pixelsize=(20E-3,20E-3))

        filenameout = os.path.join(savefolder,os.path.basename(mixedimage).split(extension)[0] + "_STED2species_MixedIntensity.tiff")
        print(filenameout)
        imsave(file=filenameout, data=imsum.astype(numpy.uint16), luts="gray", pixelsize=(20E-3,20E-3))
        filenameout = os.path.join(savefolder,
                                    os.path.basename(mixedimage).split(extension)[0] + "_STED2species_lifetimergb.tiff")

        tifffile.imwrite(filenameout, lifetime_rgb.astype(numpy.float32))

        filenameout =  os.path.join(savefolder,os.path.basename(mixedimage).split(extension)[0] + "_STED2species_f1f2.tiff")
        imagecomp=numpy.dstack((fraction2,fraction1))
        imagecomp=numpy.moveaxis(imagecomp,2,0)
        tifffile.imwrite(filenameout, imagecomp)

        fig_im3.savefig(os.path.join(savefolder, 'Images_SeparateSTED_2species_' +
                                        os.path.basename(mixedimage).split(extension)[0] + '.pdf'), transparent='True',
                        bbox_inches="tight")
    fig4.savefig(os.path.join(savefolder, 'Phasor_SeparateSTED_2species.pdf'), transparent='True',
                    bbox_inches="tight")
    fig4.savefig(os.path.join(savefolder, 'Phasor_SeparateSTED_2species_controlsonly.png'),
                    transparent='True', bbox_inches="tight")
#plt.show()

