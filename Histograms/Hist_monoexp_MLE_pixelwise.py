

""""    Computes lifetime with MLE method of each pixel
        in selected region
        Shows lifetime image modulated in intensity
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
dossier = os.path.expanduser("~/Documents/Github/2-Species-SPLIT-STED/Functions")
path1.append(dossier)
from Mono_fit import ExpFunMono_MLE
from Main_functions import load_image, get_foreground, select_channel
from tiffwrapper import LifetimeOverlayer
# -----------------------------------------------------------

filename =easygui.diropenbox(default=os.path.expanduser("~Desktop"))
print(filename)

mapcomp = {'CONF561': 'Confocal_561 {11}', 'STED561' : 'STED 561 {11}'}
mapcomp = {'CONF561': 0, 'STED561' : 1}

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

# -----------------------------------------------------------
#     Choose msr file in folder

for imagei in images:
    print(os.path.basename(imagei))
# Ask user to choose an image file
numim = int(input('Fichier msr a extraire (1er=0): '))
images = images[numim]
# Read the selected image file
imagemsr = load_image(images)



# -----------------------------------------------------------
#     Open mapcomp's images

for key in mapcomp:
    print(mapcomp[key])
    image1=select_channel(imagemsr, mapcomp[key])
    #image1=imagemsr[mapcomp[key]]
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
    indice=20
    mlifetime = numpy.empty((imsum.shape))

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

    fig1, ax1 = plt.subplots()
    imgplot1 = ax1.imshow(mlifetime, cmap='jet', vmin=0, vmax=4)
    ax1.axis('off')
    cbar =fig1.colorbar(imgplot1)
    cbar.set_label("Lifetime [ns]")
    #plt.show()

# -----------------------------------------------------------
#    Lifetime modulated in intensity


    overlayer = LifetimeOverlayer(mlifetime, imsum/imsum.max(), cname='jet')
    lifetime_rgb, cmap = overlayer.get_overlay(
        lifetime_minmax=(0, 4),
        intensity_minmax=(0, 0.7)
                )
            
    """ If error, watch out conversion in lifetime.py """

    fig2, ax2 = plt.subplots()
    ax2.axis('off')
    img = ax2.imshow(lifetime_rgb)
    cbar = fig2.colorbar(cmap, ax=ax2)
    cbar.set_label("temps de vie [ns]")
    filenameout = os.path.basename(images).split(extension)[0] + "_MLE_Lifetime_{}.tiff".format(key)
    tifffile.imwrite(filenameout, mlifetime.astype(numpy.uint16))
    filenameout = os.path.basename(images).split(extension)[0] + "_Intensity_{}.tiff".format(key)
    tifffile.imwrite(filenameout, imsum.astype(numpy.uint16))
    fig2.savefig(os.path.basename(images).split(extension)[0] +'MLELifetime_IntensityComposite_{}.png'.format(key), transparent='True', bbox_inches="tight")
    fig2.savefig(os.path.basename(images).split(extension)[0] +'MLELifetime_IntensityComposite_{}.pdf'.format(key), transparent='True', bbox_inches="tight")

    plt.show()



