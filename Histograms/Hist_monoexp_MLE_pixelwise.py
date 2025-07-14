

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
from Mono_fit import ExpFunMono_MLE
from Main_functions import load_image, get_foreground, select_channel
from tiffwrapper import LifetimeOverlayer
# -----------------------------------------------------------

filename =easygui.diropenbox(default=os.path.expanduser("~Desktop"),title="Select the folder containing the images")
print(filename)

mapcomp = {'CONF561': 'Confocal_561 {11}', 'STED561' : 'STED 561 {11}'}
#mapcomp = {'CONF561': 0, 'STED561' : 1} # For Tiff file, give the channel numbers

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

# -----------------------------------------------------------
# Ask user to choose an image file

for i,imagei in enumerate(images):
    print(i, os.path.basename(imagei))
numim = int(input('Enter the index of the file to load (1st=0): '))
images = images[numim]
# Read the selected image file
imagemsr = load_image(images)

# -----------------------------------------------------------

for key in mapcomp:
    print(mapcomp[key])
    image1=select_channel(imagemsr, mapcomp[key])
    #image1=imagemsr[mapcomp[key]]
    dim = image1.shape
    print(dim)
    centerx=int(dim[0]/2)
    centery = int(dim[1] / 2)

    if dim[2] > 250 :
        image1 = image1[:,:,:250].astype(numpy.int16)
        print(image1.shape)
        dim=image1.shape
    dim = image1.shape
    imsum= numpy.sum(image1[:,:,10:], axis=2, dtype = numpy.int16)

# -----------------------------------------------------------

    seuil= 5
    indice=20
    mlifetime = numpy.empty((imsum.shape))


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

    fig1, ax1 = plt.subplots()
    imgplot1 = ax1.imshow(mlifetime, cmap='jet', vmin=0, vmax=4)
    ax1.axis('off')
    cbar =fig1.colorbar(imgplot1)
    cbar.set_label("Lifetime [ns]")


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

    # Save the figures and images
    filenameout = os.path.basename(images).split(extension)[0] + "_MLE_Lifetime_{}.tiff".format(key)
    tifffile.imwrite(os.path.join(savefolder,filenameout), mlifetime.astype(numpy.uint16))
    filenameout = os.path.basename(images).split(extension)[0] + "_Intensity_{}.tiff".format(key)
    tifffile.imwrite(os.path.join(savefolder,filenameout), imsum.astype(numpy.uint16))
    fig1.savefig(os.path.join(savefolder,os.path.basename(imagei).split(extension)[0] +'MLELifetime_{}.pdf'.format(key)), transparent='True', bbox_inches="tight")
    fig2.savefig(os.path.join(savefolder,os.path.basename(images).split(extension)[0] +'MLELifetime_IntensityComposite_{}.png'.format(key)), transparent='True', bbox_inches="tight")
    fig2.savefig(os.path.join(savefolder,os.path.basename(images).split(extension)[0] +'MLELifetime_IntensityComposite_{}.pdf'.format(key)), transparent='True', bbox_inches="tight")

    plt.show()



