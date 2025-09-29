

""""  
       Takes an intensity image and a mean lifetime map image and returns the lifetime image modulated in intensity
"""
import os
import matplotlib.pyplot as plt
import glob
import numpy
import easygui

import tifffile
import os.path
from sys import path as path1; 
Functionspath=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Functions")
path1.append(Functionspath)

from Main_functions import load_image, get_foreground, select_channel
from tiffwrapper import LifetimeOverlayer
# -----------------------------------------------------------
# Path to the folder containing the images
# Opens file dialog to select folder
images =easygui.fileopenbox(default=os.path.expanduser("~Desktop"),title="Select the files containing the images",multiple=True)
print(images)

seuil=3

# Ask the user for a name for the output folder and create it
savefolder=str(input("Name of Output folder: "))
savefolder = os.path.join(os.path.expanduser("~/Desktop"), savefolder)
os.makedirs(savefolder, exist_ok=True)

# -----------------------------------------------------------
# Ask user to choose an image file

for i,imagei in enumerate(images):
    print(i, os.path.basename(imagei))
numim = int(input('Enter the index of the Intensity Image (1st=0): '))

# Read the selected image file
intensity_image = tifffile.imread(images[numim])
numim = int(input('Enter the index of the Lifetime Image (1st=0): '))
lifetime_image = tifffile.imread(images[numim])

lifetime_nan=lifetime_image.copy()
print('Intensity min, max: ', intensity_image.min(), intensity_image.max(),intensity_image.dtype)
print('Lifetime min, max: ', lifetime_nan.min(), lifetime_nan.max(),lifetime_nan.dtype)
lifetime_nan[intensity_image < seuil]=5000
lifetime_nan[lifetime_nan== 5000]=numpy.nan

# Make a plot of the lifetime image
fig1, ax1 = plt.subplots()
imgplot1 = ax1.imshow(lifetime_nan, cmap='jet',vmin=1, vmax=4)
ax1.axis('off')
cbar =fig1.colorbar(imgplot1)
cbar.set_label("Lifetime [ns]")

# -----------------------------------------------------------
#    Lifetime modulated in intensity


overlayer = LifetimeOverlayer(lifetime_image, intensity_image/intensity_image.max(), cname='jet')
lifetime_rgb, cmap = overlayer.get_overlay(
    lifetime_minmax=(1, 4),
    intensity_minmax=(0, 0.3)
            )
            
""" If error, watch out conversion in lifetime.py """

fig2, ax2 = plt.subplots()
ax2.axis('off')
img = ax2.imshow(lifetime_rgb)
cbar = fig2.colorbar(cmap, ax=ax2)
cbar.set_label("temps de vie [ns]")

    # Save the figures and images

fig2.savefig(os.path.join(savefolder,'IntensityComposite.pdf'), transparent='True', bbox_inches="tight")
fig1.savefig(os.path.join(savefolder,'LifetimeImage.pdf'), transparent='True', bbox_inches="tight")
plt.show()



