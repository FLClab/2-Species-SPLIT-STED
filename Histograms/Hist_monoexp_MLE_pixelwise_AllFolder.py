

""""    Computes lifetime with by performing a fit of a mono-exponential function using the MLE method for each pixel
   
       Returns the lifetime image modulated in intensity
"""

import os
import matplotlib.pyplot as plt
import glob
import numpy
import easygui
from scipy.optimize import minimize
from tqdm import tqdm
import tifffile
import os.path
from sys import path as path1; 
Functionspath=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Functions")
path1.append(Functionspath)
from Main_functions import get_foreground,load_msr,load_image,select_channel
from Mono_fit import ExpFunMono_MLE
from tiffwrapper import LifetimeOverlayer


# -----------------------------------------------------------

# Path to the folder containing the images

filename =easygui.diropenbox(default=os.path.expanduser("~Desktop"),title="Select the folder containing the images")
print(filename)

# Dictionary of the image identifiers (Channel names) to be included

mapcomp = {'STED635': 'STED_635P {2}','Conf635': 'Conf_635P {2}'}
mapcomp = {'CONF561': 'Confocal_561 {11}', 'STED561' : 'STED 561 {11}'}
#mapcomp = {'CONF': 0, 'STED' : 1} # For Tiff file, give the channel numbers

# Make list of all the images in the folder
extension = ".msr"
path=os.path.join(filename,"*.msr")
images = glob.glob(path)
print('There are ',len(images), ' msr files in this folder')
if len(images) == 0:
    path=os.path.join(filename,"*.tiff")
    images = glob.glob(path)
    print('There are ',len(images), ' tiff files in this folder')
    extension = ".tiff"

# Ask the user for a name for the output folder and create it
savefolder=str(input("Name of Output folder: "))
savefolder = os.path.join(os.path.expanduser("~/Desktop"), savefolder)
os.makedirs(savefolder, exist_ok=True)
#plt.style.use('dark_background')

# -----------------------------------------------------------
# Loop over the images in the folder
for imagei in images:
    print(os.path.basename(imagei))
    imagemsr = load_image(imagei)
    
# -----------------------------------------------------------
#  # Loop over the channels in the dictionary

    for key in mapcomp:
        print(mapcomp[key])
        #image1=imagemsr[mapcomp[key]]
        image1 = select_channel(imagemsr, mapcomp[key])
        dim = image1.shape
        print(dim)
        centerx=int(dim[0]/2)
        centery = int(dim[1] / 2)
        #image1=image1[centerx-30:centerx+30,centery-30:centery+30,:]
        if dim[2] > 250 :
            image1 = image1[:,:,:250].astype(numpy.int16)
            print(image1.shape)
            dim=image1.shape

        #image1 = image1[  140:240 , 240:340, :250 ] # Change image dimensions here
        dim = image1.shape

        imsum= numpy.sum(image1[:,:,10:], axis=2, dtype = numpy.int16)


    # -----------------------------------------------------------

        seuil= 5
        mlifetime = numpy.empty((imsum.shape))
        indice=20
        # Loop over the pixels in the image
        # For each pixel, calculate the mono-exponential fit using MLE
        # If the pixel intensity is below the threshold, set the lifetime to 0
        mlifetime.fill(numpy.nan)
        for iy, ix in tqdm(numpy.ndindex(imsum.shape)):
            y = image1[iy, ix]
            if y.sum() < seuil :
                mlifetime[iy, ix] = 0
            else :
                #maxy = numpy.max(y) 
                #indice = numpy.argmax(y)
                y = y[indice:]
                y= y / y.sum()
                absci = numpy.linspace(0,y.shape[0]-1, num =y.shape[0])*0.08
    #   Computes lifetime for each pixel

                tau = 2
                nb_photons = 100
                bounds = [(0,0.0001),(0.1,5),(0,900000)] # (min, max): Background, tau, amp
                w = minimize(ExpFunMono_MLE, 
                            x0 = [0,tau,nb_photons], 
                            args = [absci, y], 
                            bounds = bounds)
                bg, tau, amp =  w["x"]
                mlifetime[iy, ix] = w['x'][1]

        mlifetime_nan=mlifetime.copy()
        mlifetime_nan[imsum < seuil]=numpy.nan

        # Make a plot of the lifetime image
        fig1, ax1 = plt.subplots()
        imgplot1 = ax1.imshow(mlifetime_nan, cmap='jet',vmin=1, vmax=4)
        ax1.axis('off')
        cbar =fig1.colorbar(imgplot1)
        cbar.set_label("Lifetime [ns]")
        

    # -----------------------------------------------------------
    #    Lifetime modulated in intensity
        overlayer = LifetimeOverlayer(mlifetime, imsum/imsum.max(), cname='jet')
        lifetime_rgb, cmap = overlayer.get_overlay(
            lifetime_minmax=(1,4),
            intensity_minmax=(0, 0.45)
                    )
                
        # Make a plot of the lifetime image modulated in intensity

        fig2, ax2 = plt.subplots()
        ax2.axis('off')
        img = ax2.imshow(lifetime_rgb)
        cbar = fig2.colorbar(cmap, ax=ax2)
        cbar.set_label("Lifetime [ns]")

        # Save the figures and images
        filenameout = os.path.join(savefolder,os.path.basename(imagei).split(extension)[0] + "_MLE_Lifetime_{}.tiff".format(key))
        tifffile.imwrite(filenameout, mlifetime)
        filenameout =os.path.join(savefolder, os.path.basename(imagei).split(extension)[0] + "_Intensity_{}.tiff".format(key))
        tifffile.imwrite(filenameout, imsum.astype(numpy.uint16))
        fig1.savefig(os.path.join(savefolder,os.path.basename(imagei).split(extension)[0] +'MLELifetime_{}.pdf'.format(key)), transparent='True', bbox_inches="tight")
        fig2.savefig(os.path.join(savefolder,os.path.basename(imagei).split(extension)[0] +'MLELifetime_IntensityComposite_{}.pdf'.format(key)), transparent='True', bbox_inches="tight")
        plt.close('all')




