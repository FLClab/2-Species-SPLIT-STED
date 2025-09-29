"""

Program that adds synthetic background to single species STED-FLIM images and calculates the reference triangle for 2 species-SPLIT-STED.
Synthetic background is added by sampling pixel positions from a Uniform distribution and placing them in time bins by sampling from different temporal distributions.
The script then saves the triangle vertex coordinates for each noise condition in a CSV file.

The routine is defined as a function and called at the end of the script in a loop over the different STED powers and Noise quantities.

"""


import os
from random import random
from sys import path as path1


Functionspath=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Functions")
path1.append(Functionspath)
from Main_functions import (load_image,select_channel,line_equation, to_polar_coord, polar_to_cart, get_foreground)
from Phasor_functions import Median_Phasor,DTCWT_Phasor,unmix3species
from objectives import (Squirrel, Bleach)
import decorr

from tiffwrapper import imsave,LifetimeOverlayer
import math
import random
import numpy
import glob
import itertools
import tifffile
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn
import pandas as pd
import easygui
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression
import skimage

from skimage import filters
import scipy
import shapely
from skspatial.objects import Circle
from skspatial.objects import Line
#plt.style.use('dark_background')
matplotlib.rcParams['axes.linewidth'] = 0.8
# ------------------ Default Input variables ----------------
params_dict = {
    # Parameter in option in the matlab code
    "smooth_factor": 0.2,  # % 'Smoothing factor:'
    "phasor_smooth_cycles": 1,  # % 'Smoothing cycles phasor:'
    "foreground_threshold": 10,
    "harm1": 1,
}

Noiselist=[[0,0,"Uniform"]]
types=["Uniform","IRF","Alexa647"]
pixels=[1,5,10,15,20,25,50]
photons=[10,20,30,40,50]
combos=list(itertools.product(pixels,photons))
for type in types:
    print(type)
    for combo in combos:
        
        combo=list(combo)
        print("combo",combo)
        combo.append(type)
        print("combo",combo)
        Noiselist.append(combo)
print(Noiselist)
#,[1,10,"Uniform"],[1,20,"Uniform"],[1,50,"Uniform"],[1,70,"Uniform"],[5,10,"Uniform"],[10,10,"Uniform"],[15,10,"Uniform"],[20,10,"Uniform"],[25,10,"Uniform"],[1,10,"IRF"],[1,20,"IRF"],[1,50,"IRF"],[5,10,"IRF"],[10,10,"IRF"],[15,10,"IRF"],[20,10,"IRF"],[1,10,"Alexa647"],[1,20,"Alexa647"],[1,50,"Alexa647"],[5,10,"Alexa647"],[10,10,"Alexa647"],[15,10,"Alexa647"],[20,10,"Alexa647"]]# List of noise levels to simulate [percent pixels, number of photons]
#Opens dialog box for the user to select the folders containing the control images
filename1=easygui.diropenbox(default=os.path.expanduser("~Desktop"),title="Select folder containing control images for First fluorophore")
filename2=easygui.diropenbox(default=os.path.expanduser("~Desktop"),title="Select folder containing control images for Second fluorophore")


labels = ['Bassoon_CF594 Confocal','Bassoon_CF594 STED 10%','Bassoon_CF594 STED 20%','Bassoon_CF594 STED 30%','PSD95_STORANGE Confocal','PSD95_STORANGE STED 10%','PSD95_STORANGE STED 20%','PSD95_STORANGE STED 30%', 'Mixture']


# Channels to use for the control images
keys = ['Confocal_561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'Confocal_561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'STED 561 {11}']
#keys= ['Conf_635P {2}','STED_635P {2}','STED_635P {2}','STED_635P {2}','STED_635P {2}','Conf_635P {2}','STED_635P {2}','STED_635P {2}','STED_635P {2}','STED_635P {2}']
#keys = [ 'Conf640 {10}', 'STED640 {10}', 'STED640 {10}', 'STED640 {10}','Conf640 {10}', 'STED640 {10}', 'STED640 {10}', 'STED640 {10}', 'STED640 {10}']
#keys=[0,1,1,1,0,1,1,1] # For Tiff file, use channel ID

# Create a folder to save the results
savefolder = "Simulation_BuildTriangle_Noise"
savefolder=os.path.join(os.path.expanduser("~/Desktop"),savefolder)
os.makedirs(savefolder,exist_ok=True)

# -----------------------------------------------------------
def BuildTriangleWithNoise(noise):
    """
    Function that builds 2-species SPLIT-STED triangle after adding noise to the control images.

    Args:
        
        noise (float): Standard deviation of the Gaussian noise to add to the control images

    Returns:
        An array containing the coordinates of the triangle vertices for the noise level.

    Saves : 
    For each pair of images:
        - A pdf file of the phasor plot colorcoded by the mixture fraction
        - ground truth and predicted fraction images
    Overall:
        - A CSV file of the coordinates of the triangle vertices for each noise level
        - A legend file containing the paths to the control images used
    """
    
    f2= os.path.join(filename2,"*.msr")
    f1= os.path.join(filename1,"*.msr")
    filenamescontrol = [f1,f1,f1,f1, f2,f2,f2,f2]


    NUMIM=[7,7,1,0,0,0,1,19]
        # Create a dataframe to store the results
    if noise[2]=="Uniform":
        distribution=numpy.ones((1,250),dtype=int) # Uniform distribution of photons across time bins
        distribution=distribution[0,:]

    else:
        if noise[2]=="IRF":
            
            filename=r'C:\Users\FLCLab\Documents\GitHub\2-Species-SPLIT-STED\Acquisition\IRF_Measurement\GoldBead_individuel_CONF_640_Cy5_45us_50p.msr'
            key= 'STAR 635P_CONF {0}'
            #key=0 # For Tiff file,
        elif noise[2]=="Alexa647":
            filename=r'C:\Users\FLCLab\Documents\GitHub\2-Species-SPLIT-STED\Simulation\alphaTubulin_AF647_STEDPowerBleach_5to20_1_18_5percentSTED.tiff'
            key= 0 # For Tiff file,

        imagemsr = load_image(filename)
        image1=select_channel(imagemsr, key)
        dim = image1.shape
        imsum= numpy.sum(image1, axis=2)
        seuil = get_foreground(imsum)
        y = image1[imsum > seuil, :] 
        distribution = numpy.sum(y, axis=0)
        print(distribution.shape)
            


    colors=['lightsteelblue', 'deepskyblue', 'royalblue','midnightblue','lightsalmon','lightcoral','crimson','darkred','springgreen']



# Use the control images to build the triangle in phasor space that will be used to unmix the mixed images
    msrfiles = []
    for k,filename in enumerate(filenamescontrol) :
        print(labels[k])
        # Make list of all the images in the folder
        images = glob.glob(filename)
        print("Found {} images in {}".format(len(images), filename))
        # Select the image to use as control
        numim=NUMIM[k]
        image = images[numim]
        msrfiles.append(image)
    print(msrfiles)


    CoM_x, CoM_y = [], []
    # Create a figure for the phasor plot
    fig4,ax_scatter = plt.subplots(figsize=(3,3))
    # draw universal semi-circle
    edge = numpy.linspace(start=0, stop=15, num=200)
    theta = numpy.linspace(0, numpy.pi, 100)
    r = 0.5
    x1 = r * numpy.cos(theta) + 0.5
    x2 = r * numpy.sin(theta)
    ax_scatter.plot(x1, x2, color="k", ls="--",linewidth=0.8)

# Create legend file 
    with open(os.path.join(savefolder,'legend.txt'),'w') as data: 
        data.write("Controls\n")
    scatterlist=[]

    for i, msr in enumerate(msrfiles) : 
        df = pd.DataFrame(columns=['x','y'])
        dg = pd.DataFrame(columns=['g', 's'])

       # Write the control images info in the legend file
        with open(os.path.join(savefolder,'legend.txt'),'a') as data: 
            data.write("{}\t{}\t{}\n".format(labels[i],keys[i],msr))
         # Load the image and select the channel    
        imagemsr=load_image(msr)
        image1 = select_channel(imagemsr, keys[i])

        # Add Noise to the image
        print("Adding {} photons to {}percent of the pixels".format(noise[1], noise[0]))
        values=numpy.linspace(0,image1.shape[2],num=image1.shape[2],endpoint=False)
        print(values.shape)
        print(distribution.shape)
        num_noisypixels= int(noise[0]/100*image1.shape[0]*image1.shape[1]) # Calculate number of pixels to add noise to
        x = numpy.random.randint(0,high= image1.shape[0],size=num_noisypixels)
        y = numpy.random.randint(0, high=image1.shape[1],size=num_noisypixels)
       
        for n in range(num_noisypixels):
            #tlist= numpy.random.randint(0,high=image1.shape[2],size=noise[1]) # Randomly select time bins to add photons to
            

            tlist=random.choices(values, weights=distribution, k=noise[1])
            for t in tlist:
                t=int(t)
                image1[x[n], y[n],t] += 1

        imsum = image1.sum(axis=2)
        imsum = imsum.astype('int16')
        
        # Calculate the phasor distribution of the foreground of the control image
        seuil = get_foreground(imsum)
        print("Caclulation for an image of shape", image1.shape, "...")
        params_dict["foreground_threshold"] = seuil
        
        print("foreground_threshold=", params_dict["foreground_threshold"])
        
        x,y,g_smoothed,s_smoothed, original_idxes= Median_Phasor(image1, params_dict, **params_dict)
        df['x']=x.flatten()
        df['y']=y.flatten()
        m, phi = to_polar_coord(df['x'], df['y'])
        g,s =polar_to_cart(m, phi)
        dg['g'], dg['s'] = g, s

        # Find the centroid of the phasor distribution using KMeans clustering
        kmeans = KMeans(n_clusters = 1, init = 'k-means++', random_state = 42)
        y_kmeans = kmeans.fit_predict(dg)
        CoM_x.extend(kmeans.cluster_centers_[:, 0][:].tolist())
        CoM_y.extend(kmeans.cluster_centers_[:, 1][:].tolist())

        # Plot the phasor distribution
        #a=ax_scatter.scatter(g, s, s=1, c=colors[i], alpha=0.10,label=labels[i],rasterized=True)
        #scatterlist.append(a)

    ax_scatter.set_xlim(0, 1)
    ax_scatter.set_ylim(0, 1)

    ax_scatter.set_xlabel('g')
    ax_scatter.set_ylabel('s')

    ##Calculating the points of the triangle in phasor space to be used for unmixing

    # Projecting the centroids of the controls on the semi-circle
    xaxis = numpy.linspace(0, 1.5, 100)
    r = 0.5
    norm = numpy.sqrt((CoM_x[0] - 0.5) ** 2 + (CoM_y[0] ** 2))
    Pn_x = 0.5 + (r * (CoM_x[0] - 0.5) / norm)
    Pn_y = 0 + r * (CoM_y[0] - 0) / norm
    P_n = numpy.array([Pn_x, Pn_y])
    norm = numpy.sqrt((CoM_x[4] - 0.5) ** 2 + (CoM_y[4] ** 2))
    P2_x = 0.5 + (r * (CoM_x[4] - 0.5) / norm)
    P2_y = 0 + r * (CoM_y[4] - 0) / norm
    p2 = numpy.array([P2_x, P2_y])

    #Fit linear trajectory through centroids of controls for each species
    PointsSpecies1 = numpy.stack(
        [numpy.array([Pn_x, Pn_y]), numpy.array([CoM_x[1], CoM_y[1]]), numpy.array([CoM_x[2], CoM_y[2]]),
            numpy.array([CoM_x[3], CoM_y[3]])])
    PointsSpecies2 = numpy.stack(
        [numpy.array([P2_x, P2_y]), numpy.array([CoM_x[5], CoM_y[5]]), numpy.array([CoM_x[6], CoM_y[6]]),
            numpy.array([CoM_x[7], CoM_y[7]])])

    coeffs1 = numpy.polyfit([Pn_x, CoM_x[1], CoM_x[2], CoM_x[3]], [Pn_y, CoM_y[1], CoM_y[2], CoM_y[3]], 1)
    coeffs2 = numpy.polyfit([P2_x, CoM_x[5], CoM_x[6], CoM_x[7]], [P2_y, CoM_y[5], CoM_y[6], CoM_y[7]], 1)

    #Find intersection point of the two lines
    y1 = coeffs1[0] * xaxis + coeffs1[1]
    y2 = coeffs2[0] * xaxis + coeffs2[1]
    det = coeffs2[0] - coeffs1[0]
    x = (coeffs1[1] - coeffs2[1]) / det
    y = (coeffs2[0] * coeffs1[1] - coeffs1[0] * coeffs2[1]) / det
    p0 = numpy.array([x, y])
    print('p0', p0)

# Check if intersection point is inside the semi-circle
    circ = Circle((0.5, 0), radius=0.5)
    check = circ.contains_point([x, y])
    if check == False:  # If intersection point is outside the circle, find intersection with circle
        print("I'm outside the circle, coming in!")
        circle = Circle([0.5, 0], 0.5)
        line1 = Line.from_points([Pn_x, Pn_y], [x, y])
        line2 = Line.from_points([P2_x, P2_y], [x, y])
        point_a, point_b = circle.intersect_line(line1)
        point_c, point_d = circle.intersect_line(line2)

        print(numpy.array([point_a, point_c]))
        p0 = numpy.mean(numpy.array([point_a, point_c]), axis=0)
        print("New p0", p0)
        x = p0[0]
        y = p0[1]
    if y < 0:  # If intersection point is under the semi-circle, find intersection with x-axis
        print("I'm in the negatives, coming up!")
        line1 = Line.from_points([Pn_x, Pn_y], [x, y])
        line2 = Line.from_points([P2_x, P2_y], [x, y])
        line3 = Line.from_points([0, 0], [1,0])
        point_a = line3.intersect_line(line1)
        point_b = line3.intersect_line(line2)
        print(numpy.array([point_a, point_b]))
        p0 = numpy.mean(numpy.array([point_a, point_b]), axis=0)
        print("New p0", p0)
        x = p0[0]
        y = p0[1]

    print('p0', p0)
    print("POINTS", P_n, p2, p0)

# Plot the lines and points in the phasor graph
    centroidscatter1=ax_scatter.scatter(PointsSpecies1[:,0],PointsSpecies1[:,1],s=50,c="firebrick")
    centroidscatter2 = ax_scatter.scatter(PointsSpecies2[:,0],PointsSpecies2[:,1], s=50, c="firebrick")
    pnscatter = ax_scatter.scatter(Pn_x, Pn_y, s=50, c='gold')
    p2scatter = ax_scatter.scatter(P2_x, P2_y, s=50, c='gold')
    p0scatter = ax_scatter.scatter(p0[0], p0[1], s=50, c='gold')
    p0pnline = ax_scatter.plot([Pn_x, p0[0]], [Pn_y, p0[1]], c='dodgerblue')
    p2pnline = ax_scatter.plot([Pn_x, P2_x], [Pn_y, P2_y], c='dodgerblue')
    p0p2line = ax_scatter.plot([P2_x, p0[0]], [P2_y, p0[1]], c='dodgerblue')

    

    # Save the phasor plot
    fig4.savefig(os.path.join(savefolder, "Phasor_ControlsOnly_{}Noise_{}photons_{}percentpixels.pdf".format(noise[2],noise[1], noise[0])), transparent='True',
                 bbox_inches="tight",dpi=900)
    
    # Remove the previous scatter plots from the figure
    t = [scatter.remove() for scatter in scatterlist]
    pnscatter.remove()
    p2scatter.remove()
    p0scatter.remove()
    centroidscatter1.remove()
    centroidscatter2.remove()
    ax_scatter.lines[-1].remove()
    ax_scatter.lines[-1].remove()
    ax_scatter.lines[-1].remove()




    return [noise[0],noise[1],noise[2],Pn_x, Pn_y, P2_x, P2_y, p0[0], p0[1]]


# Main code to run the simulation
Overall_data=pd.DataFrame(columns=["Noise percent pixels","Noise number photons","Noise type",'P1_x', 'P1_y', 'P2_x', 'P2_y', "P3_x", "P3_y"])



for n, noise in enumerate(Noiselist):
    stats=BuildTriangleWithNoise(noise)
    Overall_data.loc[n] = stats



for index, row in Overall_data.iterrows():
    if row["Noise percent pixels"]==0 and row["Noise number photons"]==0:
        print("I am control")
        control = shapely.Polygon(((row['P1_x'], row['P1_y']), (row['P2_x'], row['P2_y']), (row['P3_x'], row['P3_y'])))
    polygon = shapely.Polygon(((row['P1_x'], row['P1_y']), (row['P2_x'], row['P2_y']), (row['P3_x'], row['P3_y'])))
    i = control.intersection(polygon)
    u=control.union(polygon)
    iou=i.area/u.area
    print("IoU for noise level {}% pixels, {} photons: {}".format(row["Noise percent pixels"], row["Noise number photons"], iou))
    Overall_data.loc[index,"IoU"]=iou
Overall_data.to_csv(os.path.join(savefolder,"Overall_data_3species_Noise.csv"), index=False)







