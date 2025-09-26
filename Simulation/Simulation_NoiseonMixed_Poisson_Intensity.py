"""

Program that creates synthetic 2 species STED-FLIM images by summing single species images and then 
unmixing them using the 2 Species SPLIT-STED method.
The program then compares the unmixed images to the ground truth images and computes different metrics such as resolution and nanoJ-SQUIRREL


The routine is defined as a function and called at the end of the script in a loop over the different STED powers.

"""


import os
from sys import path as path1;
Functionspath=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Functions")
path1.append(Functionspath)
from Main_functions import (load_image,select_channel,line_equation, to_polar_coord, polar_to_cart, get_foreground)
from Phasor_functions import Median_Phasor,DTCWT_Phasor,unmix3species
from objectives import (Squirrel, Bleach)
import decorr

from tiffwrapper import imsave,LifetimeOverlayer
import math

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
import random
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


def abc_to_rgb(A=0.0,B=0.0,C=0.0):
# Map values A, B, C (all in domain [0,1]) to
# suitable red, green, blue values.
    return (min(B+C,1.0),min(A+C,1.0),min(A+B,1.0))

#Opens dialog box for the user to select the folders containing the control images
filename1=easygui.diropenbox(default=os.path.expanduser("~Desktop"),title="Select folder containing control images for First fluorophore")
filename2=easygui.diropenbox(default=os.path.expanduser("~Desktop"),title="Select folder containing control images for Second fluorophore")


# Image IDs to use as controls for each STED power
#Powerslist=[[10,[7,7,1,0,0,0,1,19]],[20,[7,7,1,0,0,0,1,19]],[30,[7,7,1,0,0,0,1,19]],[40,[7,7,1,0,0,0,1,19]]] #PSD95 Bassoon Cy3
Powerslist=[[20,[7,7,1,0,0,0,1,19]],[40,[7,7,1,0,0,0,1,19]]] #PSD95 Bassoon Cy3
Powerslist=[[10,[3,3,1,0,0,0,1,8]],[20,[3,3,1,0,0,0,1,8]],[30,[3,3,1,0,0,0,1,8]],[40,[3,3,1,0,0,0,1,8]]] #PSD95 Bassoon Cy3 MiniNew_10,20,30
#Powerslist=[[20,[2,2,1,0,0,0,1,4]]] #PSD95 Bassoon Cy3 in MiniNew
#Powerslist=[[5,[0,0,1,9,22,22,7,0]],[10,[0,0,1,9,22,22,7,0]],[15,[0,0,1,9,22,22,7,0]],[20,[0,0,1,9,22,22,7,0]]]# Spectrin Bassoon Cy5

labels = ['Bassoon_CF594 Confocal','Bassoon_CF594 STED 10%','Bassoon_CF594 STED 20%','Bassoon_CF594 STED 30%','PSD95_STORANGE Confocal','PSD95_STORANGE STED 10%','PSD95_STORANGE STED 20%','PSD95_STORANGE STED 30%', 'Mixture']
Noiselist=[]
types=["Uniform","IRF","Alexa647"]
#pixels=[1,5,10,15,20,25,50]
#photons=[10,20,30,40,50]
pixels=[100]
photons=[3,5,7]
random.seed(42)
numpy.random.seed(seed=42)


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
Noiselist.append([0,0,"Uniform"])
print(Noiselist)
# Channels to use for the control images
keys = ['Confocal_561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'Confocal_561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'STED 561 {11}']
#keys= ['Conf_635P {2}','STED_635P {2}','STED_635P {2}','STED_635P {2}','STED_635P {2}','Conf_635P {2}','STED_635P {2}','STED_635P {2}','STED_635P {2}','STED_635P {2}']
#keys = [ 'Conf640 {10}', 'STED640 {10}', 'STED640 {10}', 'STED640 {10}','Conf640 {10}', 'STED640 {10}', 'STED640 {10}', 'STED640 {10}', 'STED640 {10}']
#keys=[0,1,1,1,0,1,1,1] # For Tiff file, use channel ID

# Channel to use for the mixed images
keysmixed = ['STED 561 {11}','STED 561 {11}']
#keysmixed = ['STED_635P {2}','STED_635P {2}']
#keysmixed = ['STED640 {10}', 'STED640 {10}']
#keysmixed=[1,1] # For Tiff file, use channel ID
# -----------------------------------------------------------
def Simulate3speciesLineControls(STEDPOWER, NUMIM, Noiselist):
    """
    Generates simulated combinations of pairs of images acquired with the same STED power
      and performs 2species SPLIT-STED to unmix them using control images designated by numim.

    Args:
        STEDPOWER (int): The STED power.
        NUMIM (list): A list of integers representing the control images to use for each species.

    Returns:
        An array containing the resolutions and nanoJ-SQUIRREL scores for all the pairs of images of the STED power.

    Saves : 
    For each pair of images:
        - A pdf file of the phasor plot colorcoded by the mixture fraction
        - ground truth and predicted fraction images
    Overall:
        - A CSV file of the measured metrics for all the pairs of images
        - A legend file containing the paths to the control images and to the files that are mixed together for each pairID
    """
    
    f2= os.path.join(filename2,"*.msr")
    f1= os.path.join(filename1,"*.msr")
    filenamescontrol = [f1,f1,f1,f1, f2,f2,f2,f2]
    # Create a list of the images acquired with the correct STED power
    f2= os.path.join(filename2,"*_{}percentSTED.msr".format(STEDPOWER))
    f1= os.path.join(filename1,"*_{}percentSTED.msr".format(STEDPOWER))
    filenames = [f1,f2]



    colors=['lightsteelblue', 'deepskyblue', 'royalblue','midnightblue','lightsalmon','lightcoral','crimson','darkred','springgreen']

    # Create a folder to save the results
    savefolder = "Simulation_Cy3_{}Percent_3Species_LineControls_PSD95Bassoon_NoisePoisson_Intensity_Seed42".format(STEDPOWER)
    savefolder=os.path.join(os.path.expanduser("~/Desktop"),savefolder)
    os.makedirs(savefolder,exist_ok=True)



# Combine all possible pairs of images, perform unmixing and calculate metrics
    images=[glob.glob(filename)for filename in filenames]
    number = [len(glob.glob(filename)) for filename in filenames]
    print('There are ',number, 'Images in these folders')
    
    # make list of all possible pairs of images
    pairs = list(itertools.product(images[0], images[1]))
    print(len(pairs))


    # Create a dataframe to store the results
    Overall_data=pd.DataFrame(columns=["Power",'image1', 'image2',"Foreground_int","Background_int","Foreground_int_post","Background_int_post","Noise percent pixels","Noise number photons","Noise type"])
    pairnoiseid=0
    # Loop over different noise levels
    for noise in Noiselist:
        if noise[2]=="Uniform":
            distribution=numpy.ones((1,250),dtype=int) # Uniform distribution of photons across time bins
            distribution=distribution[0,:]

        else:
            if noise[2]=="IRF":
                
                filename=r'C:\Users\FLCLab\Documents\GitHub\2-Species-SPLIT-STED\Acquisition\IRF_Measurement\GoldBead_individuel_CONF_640_Cy5_45us_50p.msr'
                #filename=r'C:\Users\andde\Documents\GitHub\2-Species-SPLIT-STED\Acquisition\IRF_Measurement\GoldBead_individuel_CONF_640_Cy5_45us_50p.msr'
                key= 'STAR 635P_CONF {0}'
                #key=0 # For Tiff file,
            elif noise[2]=="Alexa647":
                filename=r'C:\Users\FLCLab\Documents\GitHub\2-Species-SPLIT-STED\Simulation\alphaTubulin_AF647_STEDPowerBleach_5to20_1_18_5percentSTED.tiff'
                #filename=r'C:\Users\andde\Documents\GitHub\2-Species-SPLIT-STED\Simulation\alphaTubulin_AF647_STEDPowerBleach_5to20_1_18_5percentSTED.tiff'
                key= 0 # For Tiff file,

            imagemsr = load_image(filename)
            image1=select_channel(imagemsr, key)
            dim = image1.shape
            imsum= numpy.sum(image1, axis=2)
            seuil = get_foreground(imsum)
            y = image1[imsum > seuil, :] 
            distribution = numpy.sum(y, axis=0)
            print(distribution.shape)
        


        # Loop over all pairs of images
        for Pair_id,(a,b) in enumerate(pairs):
            ov_data = [STEDPOWER,a,b]
            msrfiles=[a,b]
            print("***********************************************************")
            print("Working on pair number ",Pair_id," out of ",len(pairs))
            print("***********************************************************")
            print(msrfiles)
            croplist=[]
            Imagelist = []
            Masklist = []
            Propslist = []
            Cropslist = []
            CropImageList = []
            ControlImagesList = []
            ComboMasklist = []

            seuils=[]
            # For each image in the pair, calculate the resolution and the foreground mask
            for i, msr in enumerate(msrfiles):

                imagemsr=load_image(msr)
                #print(imagemsr.keys())
                image1 = select_channel(imagemsr, keysmixed[i])
                #image1 = imagemsr[keysmixed[i]]
               
                Imagelist.append(image1)

                print(image1.shape)
                #seuil=3
                seuil=get_foreground(numpy.sum(image1[:,:,10:111],axis=2))
                seuils.append(seuil)

                mask=numpy.sum(image1[:,:,10:111],axis=2)>seuil
                mask=scipy.ndimage.binary_fill_holes(mask)
                Masklist.append(mask)

            # Create empty images of the correct size to store the combined image and masksS
            Combo=numpy.zeros(numpy.max([Imagelist[0].shape,Imagelist[1].shape],axis=0))

            print('Combo',Combo.shape)


            minx=0
            maxx=Imagelist[0].shape[0]
            miny=0
            maxy=Imagelist[0].shape[1]
            print(minx, miny, maxx, maxy)
            # Add the first image to the combined image 
            Combo[minx:maxx,miny:maxy ,:]+=Imagelist[0][:,:, :]

            minx = 0
            maxx = Imagelist[1].shape[0]
            miny = 0
            maxy = Imagelist[1].shape[1]
            print(minx, miny, maxx, maxy)
            # Add the second image to the combined image 
            Combo[minx:maxx, miny:maxy, :] += Imagelist[1][:, :, :]


            image1=Combo.copy()
            



            # Calculate the phasor of the foreground of the combined image 
            
            df = pd.DataFrame(columns=['x', 'y'])
            dg = pd.DataFrame(columns=['g', 's'])
            print(image1.shape)
            imsum = image1.sum(axis=2)
            imsum = imsum.astype('int16')
            seuil=min(seuils)
            ov_data.extend([numpy.mean(imsum[imsum>seuil]),numpy.mean(imsum[imsum<seuil])])


                    # Add Noise to the image
            print("Adding {} photons to {}percent of the pixels using {} Distribution".format(noise[1], noise[0],noise[2]))
            if noise[1]>0:
                values=numpy.linspace(0,image1.shape[2],num=image1.shape[2],endpoint=False)
                print(values.shape)
                print(distribution.shape)

                num_noisypixels= int(noise[0]/100*image1.shape[0]*image1.shape[1]) # Calculate number of pixels to add noise to
                noisegrid=numpy.random.poisson(lam=noise[1], size=(image1.shape[0], image1.shape[1])).astype(int)
                #noisegrid=int(noisegrid)
                print("noisegrid.shape:", noisegrid.shape)
                print("noisegrid min:", noisegrid.min())
                print("noisegrid max:", noisegrid.max())
                print("noisegrid mean:", noisegrid.mean())



                for x,y in numpy.ndindex(noisegrid.shape):
                    numphotons = noisegrid[x, y]
                    tlist=random.choices(values, weights=distribution, k=numphotons)
                    for t in tlist:
                        t=int(t)
                        image1[x, y,t] += 1
                #num_noisypixels= int(noise[0]/100*image1.shape[0]*image1.shape[1]) # Calculate number of pixels to add noise to
                #x = numpy.random.randint(0,high= image1.shape[0],size=num_noisypixels)
                #y = numpy.random.randint(0, high=image1.shape[1],size=num_noisypixels)
            
                #for n in range(num_noisypixels):
                    #tlist= numpy.random.randint(0,high=image1.shape[2],size=noise[1]) # Randomly select time bins to add photons to
                    

                #   tlist=random.choices(values, weights=distribution, k=noise[1])
                #  for t in tlist:
                #     t=int(t)
                    #    image1[x[n], y[n],t] += 1
            imsumpost = image1.sum(axis=2)
            imsumpost = imsumpost.astype('int16')
            ov_data.extend([numpy.mean(imsumpost[imsum>seuil]),numpy.mean(imsumpost[imsum<seuil]),noise[0], noise[1], noise[2]])

            

            Overall_data.loc[pairnoiseid] = ov_data
            pairnoiseid+=1
           

    Overall_data.to_csv(os.path.join(savefolder,"Overall_data_3species_{}.csv".format(STEDPOWER)))


    return 5


for power in Powerslist:
    stats=Simulate3speciesLineControls(power[0],power[1],Noiselist)
    #print(stats)
    plt.close("all")



