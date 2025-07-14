
"""


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
from sys import path as path1; 
Functionspath=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Functions")
path1.append(Functionspath)

from Main_functions import (to_polar_coord, polar_to_cart, get_foreground,load_image,select_channel)
from Phasor_functions import Median_Phasor
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
#    Paths to folders containing the images, 1 per dye.
filename=easygui.diropenbox(default=os.path.expanduser("~Desktop"),title="Select the folder containing the images")

#    List of keys to be used to extract the images from the msr files
keys=['Confocal_561 {11}']
keys=[0]

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

colors_Centroids=["k","m","k"]
colors=['gray', "hotpink"]

IRF= (0.9527011687260826, 0.4695955819269703)
MeanPositions={}
Ellipsedims={}



#Create the figure and axes for the 2D plot of centroids and ellipses 
fig_centroids,ax_centroids = plt.subplots(figsize=(4,4))
theta = np.linspace(0, np.pi, 100)
r = 0.5
x1 = r * np.cos(theta) + 0.5
x2 = r * np.sin(theta)
ax_centroids.plot(x1, x2, color="black", ls="--",linewidth=0.8)
ax_centroids.set_xlim(-0.02, 1.02)
ax_centroids.set_ylim(-0.02, 1.02)
ax_centroids.set_xlabel('g')
ax_centroids.set_ylabel('s')


#Read the image,calculate the phasor distribution and plot the scatter plot of the phasor distribution

# Read the image file
imagemsr = load_image(image)
print(os.path.basename(image))


# Calculate the phasor distribution for all pixels in the image that are above the foreground threshold
df = pd.DataFrame(columns=['x','y'])
dg = pd.DataFrame(columns=['g', 's'])
image1=select_channel(imagemsr, keys[0])

print("Caclulation for an image of shape", image1.shape, "...")
params_dict["foreground_threshold"]=10
print("foreground_threshold=", params_dict["foreground_threshold"])
x,y,g_smoothed,s_smoothed, orginal_idxs= Median_Phasor(image1, params_dict, **params_dict)
df['x']=x.flatten()
df['y']=y.flatten()
# Plot the scatter plot of the phasor distribution before calibration
ax_centroids.scatter(df['x'],df['y'],s=0.5, c=colors[0],alpha=0.1,rasterized=True)

# Calculate the centroid of the phasor distribution
kmeans = KMeans(n_clusters=1, init='k-means++', random_state=42)
y_kmeans = kmeans.fit_predict(df)
CoM_x=kmeans.cluster_centers_[:, 0][:]
CoM_y=kmeans.cluster_centers_[:, 1][:]

ax_centroids.scatter(CoM_x, CoM_y, color=colors_Centroids[0], s=25)

# Apply the calibration to the phasor distribution using the IRF measurement in polar coordinates and return the cartesian coordinates
m, phi = to_polar_coord(df['x'], df['y'])
g,s =polar_to_cart(m, phi)
dg['g'], dg['s'] = g, s
# Plot the scatter plot of the phasor distribution after calibration
ax_centroids.scatter(dg['g'], dg['s'],s=0.5, c=colors[1],alpha=0.1,rasterized=True)

# Plot the IRF calibration point
calibx,caliby=polar_to_cart([IRF[0],1],[IRF[1],0])
ax_centroids.scatter(calibx,caliby,color=colors_Centroids[1], s=25)


# Calculate the centroid of the calibrated phasor distribution
kmeans = KMeans(n_clusters=1, init='k-means++', random_state=42)
y_kmeans = kmeans.fit_predict(dg)
CoM_x=kmeans.cluster_centers_[:, 0][:]
CoM_y=kmeans.cluster_centers_[:, 1][:]

ax_centroids.scatter(CoM_x, CoM_y, color=colors_Centroids[2], s=25)


fig_centroids.savefig("Centroids_single.pdf",transparent=True, bbox_inches="tight",dpi=900)

plt.show()
