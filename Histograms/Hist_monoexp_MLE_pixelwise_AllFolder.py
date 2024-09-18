

""""    Computes lifetime with MLE method of each pixel
        in selected region
        Shows lifetime image modulated in intensity
"""

import os
import matplotlib.pyplot as plt
import glob
import numpy
import easygui
from sys import path as syspath; syspath.append(os.path.expanduser("~/Documents/Github/TCSPC/Analyse/Fit - MLE - TDE/"))
from Mono import fit
from scipy.optimize import minimize
from tqdm import tqdm
import tifffile
import os.path
from sys import path as path1; path1.append('/Users/marielafontaine/Documents/GitHub/Abberior-STED-FLIM/Functions')
dossier = os.path.expanduser("~/Documents/Github/Abberior-STED-FLIM/Functions")
path1.append(dossier)
from Main_functions import (choose_msr_file, get_foreground)
from convertmsr_bioformatsAB import MSRReader
from tiffwrapper import LifetimeOverlayer
# -----------------------------------------------------------

# Path to the folder containing the images

#filename =easygui.diropenbox(default=os.path.expanduser("~Desktop"))
filename = os.path.join('T:', os.sep, 'adeschenes', 'SimulationDataset_STEDFLIM', 'Cy3', "PSD95_STORANGE")
filename = os.path.join('T:', os.sep, 'adeschenes', 'SimulationDataset_STEDFLIM', 'Cy3', "rabBassoon_CF594")
filename  = os.path.join('T:', os.sep, 'adeschenes', 'SimulationDataset_STEDFLIM', 'Cy3', 'Homer_STORANGE',"MediumAcq","prefs")
print(filename)

# Dictionary of the image identifiers (Channel names) to be included

#mapcomp = {'STED635': 'STED_635P {2}','Conf635': 'Conf_635P {2}'}
mapcomp = {'CONF561': 'Confocal_561 {11}', 'STED561' : 'STED 561 {11}'}

path=os.path.join(filename,"*.msr")
images = glob.glob(path)
print('There are ',len(images), 'Images in this folder')

# Ask the user for a name for the output folder and create it
savefolder=str(input("Name of Output folder: "))
savefolder = os.path.join(os.path.expanduser("~/Desktop"), savefolder)
os.makedirs(savefolder, exist_ok=True)
#plt.style.use('dark_background')

# -----------------------------------------------------------
#     

with MSRReader() as msrreader:

    # Loop over the images in the folder
    for imagei in images:
        print(os.path.basename(imagei))
        imagemsr = msrreader.read(imagei)
        
# -----------------------------------------------------------
#  # Loop over the channels in the dictionary

        for key in mapcomp:
            print(mapcomp[key])
            image1=imagemsr[mapcomp[key]]
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
            # Make a plot of the intensity image
            fig, ax = plt.subplots()
            ax.axis('off')
            imgplot1 = ax.imshow(imsum, cmap='hot')
            cbar = fig.colorbar(imgplot1)
            cbar.set_label('Intensity')
            #plt.show()

        # -----------------------------------------------------------
        #    Lifetime matrix

            #seuil = get_foreground(image1)
            seuil= 5
            mlifetime = numpy.empty((imsum.shape))

            for iy, ix in tqdm(numpy.ndindex(imsum.shape)):
                y = image1[iy, ix]
                if y.sum() < seuil :
                    mlifetime[iy, ix] = 0
                else :
                    maxy = numpy.max(y) 
                    indice = numpy.argmax(y)
                    y = y[indice:]
                    y= y / y.sum()
                    absci = numpy.linspace(0,y.shape[0]-1, num =y.shape[0])*0.08
        #   Computes lifetime for each pixel

                    tau = 2
                    nb_photons = 100
                    bounds = [(0,0.0001),(0.1,5),(0,900000)] # (min, max): Background, tau, amp
                    w = minimize(fit.ExpFunMono_MLE, 
                                x0 = [0,tau,nb_photons], 
                                args = [absci, y], 
                                bounds = bounds)
                    bg, tau, amp =  w["x"]
                    mlifetime[iy, ix] = w['x'][1]

            mlifetime_nan=mlifetime.copy()
            mlifetime_nan[imsum < seuil]=numpy.nan
            fig1, ax1 = plt.subplots()
            imgplot1 = ax1.imshow(mlifetime_nan, cmap='jet',vmin=1, vmax=4)
            ax1.axis('off')
            cbar =fig1.colorbar(imgplot1)
            cbar.set_label("Lifetime [ns]")
            #plt.show()

        # -----------------------------------------------------------
        #    Lifetime modulated in intensity


            overlayer = LifetimeOverlayer(mlifetime, imsum/imsum.max(), cname='jet')
            lifetime_rgb, cmap = overlayer.get_overlay(
                lifetime_minmax=(1,4),
                intensity_minmax=(0, 0.45)
                        )
                    
            """ If error, watch out conversion in lifetime.py """

            fig2, ax2 = plt.subplots()
            ax2.axis('off')
            img = ax2.imshow(lifetime_rgb)
            cbar = fig2.colorbar(cmap, ax=ax2)
            cbar.set_label("temps de vie [ns]")
            filenameout = os.path.join(savefolder,os.path.basename(imagei).split(".msr")[0] + "_MLE_Lifetime_{}.tiff".format(key))
            tifffile.imwrite(filenameout, mlifetime)
            filenameout =os.path.join(savefolder, os.path.basename(imagei).split(".msr")[0] + "_Intensity_{}.tiff".format(key))
            tifffile.imwrite(filenameout, imsum.astype(numpy.uint16))
            fig1.savefig(os.path.join(savefolder,os.path.basename(imagei).split(".msr")[0] +'MLELifetime_{}.pdf'.format(key)), transparent='True', bbox_inches="tight")
            fig2.savefig(os.path.join(savefolder,os.path.basename(imagei).split(".msr")[0] +'MLELifetime_IntensityComposite_{}.pdf'.format(key)), transparent='True', bbox_inches="tight")
            plt.close('all')
   # plt.show()



