
"""
Script to compare the phasor distribution of an image using different filtering methods:
1. Raw phasor distribution 
2. Median filtered phasor distribution
3. DTCWT filtered phasor distribution

The script will plot the phasor distribution for each method and calculate the centroid of the phasor distribution.
It will also save the plots in a PDF file.

"""
import skimage

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn
import os
import itertools
import easygui
import math
import os.path
import time
from sys import path as path1; 
Functionspath=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Functions")
path1.append(Functionspath)
from Main_functions import (load_image,select_channel,to_polar_coord, polar_to_cart, get_foreground)
from Phasor_functions import Median_Phasor, DTCWT_Phasor
import scipy
from shapely.geometry.point import Point
from shapely import affinity
from sklearn.cluster import KMeans
matplotlib.rcParams['axes.linewidth'] = 0.8


# ------------------ Default Input variables ----------------
params_dict = {
  
    "smooth_factor": 0.2,  # % 'Smoothing factor:'
    "phasor_smooth_cycles": 1,  # % 'Smoothing cycles phasor:'
    "foreground_threshold": 10,
    "harm1": 1,
}
# -----------------------------------------------------------
# Path to folder containing the images
filename=easygui.diropenbox(default=os.path.expanduser("~Desktop"),title="Select the folder containing the images")


# Make list of all the images in the folder
extension = ".msr"
path = os.path.join(filename, '*.msr' )
images = glob.glob(path)
print('There are ',len(images), ' msr files in this folder')
if len(images) == 0:
    path = os.path.join(filename, '*.tiff' )
    images = glob.glob(path)
    print('There are ',len(images), ' tiff files in this folder')
    extension = ".tiff"


for i, imagei in enumerate(images):
    print(i, os.path.basename(imagei))
# Ask user to choose an image file
numim = int(input('Enter the index of the file to load (1st=0): '))
image = images[numim]

#    List of keys to be used to extract the images from the msr files
keys=['Confocal_561 {11}']
keys=[0]

colors_Centroids=["k","k","k"]
colors=[ "blue",  "blue", "blue"]


MeanPositions={}
Ellipsedims={}

#Create the figure and axes for the 2D plot of centroids and ellipses 
fig_centroids,ax_centroids = plt.subplots(figsize=(9,3),ncols=3)
theta = np.linspace(0, np.pi, 100)
r = 0.5
x1 = r * np.cos(theta) + 0.5
x2 = r * np.sin(theta)


# Read the image file
imagemsr = load_image(image)
print(os.path.basename(image))
#print(imagemsr.keys())

# Calculate the phasor distribution for all pixels in the image that are above the foreground threshold
df = pd.DataFrame(columns=['x','y'])
dg = pd.DataFrame(columns=['g', 's'])
image1=select_channel(imagemsr, keys[0])

imsum = np.sum(image1, axis=2)
print("Caclulation for an image of shape", image1.shape, "...")
#params_dict["foreground_threshold"] = get_foreground(image1)
params_dict["foreground_threshold"]=10


# Calculate the phasor distribution with and without using the median filter
for j,sm in enumerate([0,1]):
    ax_centroids[j].plot(x1, x2, color="black", ls="--",linewidth=0.8)
    ax_centroids[j].set_xlim(0, 1.)
    ax_centroids[j].set_ylim(0, 1.)
    ax_centroids[j].set_xlabel('g')
    ax_centroids[j].set_ylabel('s')


    params_dict[ "phasor_smooth_cycles"]=sm
    print("foreground_threshold=", params_dict["foreground_threshold"])
    x,y,g_smoothed,s_smoothed, orginal_idxs= Median_Phasor(image1, params_dict, **params_dict)
    df['x']=x.flatten()
    df['y']=y.flatten()
    # Apply the calibration to the phasor distribution using the IRF measurement in polar coordinates and return the cartesian coordinates
    m, phi = to_polar_coord(df['x'], df['y'])
    g,s =polar_to_cart(m, phi)
    dg['g'], dg['s'] = g, s
    ax_centroids[j].scatter(dg['g'], dg['s'],s=0.5, c=colors[j],alpha=1,rasterized=True)

    # Calculate the centroid of the phasor distribution
    kmeans = KMeans(n_clusters=1, init='k-means++', random_state=42)
    y_kmeans = kmeans.fit_predict(dg)
    CoM_x=kmeans.cluster_centers_[:, 0][:]
    CoM_y=kmeans.cluster_centers_[:, 1][:]

    ax_centroids[j].scatter(CoM_x, CoM_y, color=colors_Centroids[j], s=25)


ax_centroids[2].plot(x1, x2, color="black", ls="--",linewidth=0.8)
ax_centroids[2].set_xlim(0, 1.)
ax_centroids[2].set_ylim(0, 1)
ax_centroids[2].set_xlabel('g')
ax_centroids[2].set_ylabel('s')

# Calculate the phasor distribution with DTCWT filtering
df = pd.DataFrame(columns=['x','y'])
dg = pd.DataFrame(columns=['g', 's'])
start_DTCWT_Time = time.time()
x, y, original_indexs, Images, Images_Filtered = DTCWT_Phasor(image1,0, 2,50)
stop_DTCWT_Time = time.time()
DTCWT_Time = stop_DTCWT_Time - start_DTCWT_Time
print("DTCWT_Time", DTCWT_Time)
# Filter the phasor distribution to remove the background
x = x[imsum > params_dict["foreground_threshold"]]
y = y[imsum > params_dict["foreground_threshold"]]
Image_indices = np.arange(len(imsum.flatten())).reshape(imsum.shape)
Image_indices = Image_indices[imsum > params_dict["foreground_threshold"]]
#print('Image_indices', numpy.min(Image_indices), numpy.max(Image_indices))
df['x'] = x.flatten()
df['y'] = y.flatten()
# Apply phasor calibration in polar coordinates (based on IRF measurement) and return to cartesan (g,s)
m, phi = to_polar_coord(x.flatten(), y.flatten())
g, s = polar_to_cart(m, phi)
g,s = np.array(g),np.array(s)
indexes = np.where((g > 0) & (g < 1) & (s > 0) & (s < 1))
g,s = g[indexes],s[indexes]
original_idxes = Image_indices[indexes]
#original_idxes = Image_indices

dg['g'] = g
dg['s'] = s
ax_centroids[2].scatter(dg['g'], dg['s'],s=0.5, c=colors[2],alpha=1,rasterized=True)

# Calculate the centroid of the phasor distribution
kmeans = KMeans(n_clusters=1, init='k-means++', random_state=42)
y_kmeans = kmeans.fit_predict(dg)
CoM_x=kmeans.cluster_centers_[:, 0][:]
CoM_y=kmeans.cluster_centers_[:, 1][:]

ax_centroids[2].scatter(CoM_x, CoM_y, color=colors_Centroids[2], s=25)

fig_centroids.savefig("Phasors_CompareFilters.pdf",transparent=True, bbox_inches="tight",dpi=900)

plt.show()
