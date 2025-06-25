
"""
Code to make 3D graph of the phasor distribution of STED images of 2 different dyes for different depletion powers.
Each dye is represented by a different color scatter plot and ellipses are fit on the 70th percentile of the distributions.


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

# ------------------ Functions to plot ellipses on the 3D axes ----------------
def create_ellipse(center, lengths, angle=0):
    """
    create a shapely ellipse. adapted from
    https://gis.stackexchange.com/a/243462
    """
    circ = Point(center).buffer(1)
    print(circ.area)
    ell = affinity.scale(circ, lengths[0]/2, lengths[1]/2)

    ellr = affinity.rotate(ell, angle)

    return ellr
from mpl_toolkits.mplot3d import art3d

def rotation_matrix(d):
    """
    Calculates a rotation matrix given a vector d. The direction of d
    corresponds to the rotation axis. The length of d corresponds to
    the sin of the angle of rotation.

    Variant of: http://mail.scipy.org/pipermail/numpy-discussion/2009-March/040806.html
    """
    sin_angle = np.linalg.norm(d)

    if sin_angle == 0:
        return np.identity(3)

    d /= sin_angle

    eye = np.eye(3)
    ddt = np.outer(d, d)
    skew = np.array([[    0,  d[2],  -d[1]],
                  [-d[2],     0,  d[0]],
                  [d[1], -d[0],    0]], dtype=np.float64)

    M = ddt + np.sqrt(1 - sin_angle**2) * (eye - ddt) + sin_angle * skew
    return M

def pathpatch_2d_to_3d(pathpatch, z = 0, normal = 'z'):
    """
    Transforms a 2D Patch to a 3D patch using the given normal vector.

    The patch is projected into they XY plane, rotated about the origin
    and finally translated by z.
    """
    if type(normal) is str: #Translate strings to normal vectors
        index = "xyz".index(normal)
        normal = np.roll((1.0,0,0), index)

    normal /= np.linalg.norm(normal) #Make sure the vector is normalised

    path = pathpatch.get_path() #Get the path and the associated transform
    trans = pathpatch.get_patch_transform()

    path = trans.transform_path(path) #Apply the transform

    pathpatch.__class__ = art3d.PathPatch3D #Change the class
    pathpatch._code3d = path.codes #Copy the codes
    pathpatch._facecolor3d = pathpatch.get_facecolor #Get the face color

    verts = path.vertices #Get the vertices in 2D

    d = np.cross(normal, (0, 0, 1)) #Obtain the rotation vector
    M = rotation_matrix(d) #Get the rotation matrix

    pathpatch._segment3d = np.array([np.dot(M, (x, y, 0)) + (0, 0, z) for x, y in verts])

def pathpatch_translate(pathpatch, delta):
    """
    Translates the 3D pathpatch by the amount delta.
    """
    pathpatch._segment3d += delta
# ------------------ Default Input variables ----------------
params_dict = {
    # Parameter in option in the matlab code
    #    "Tg" : 6, #% 'First frame to sum:'
    "Nb_to_sum": 250,  # The Tg infered from this variable override Tg
    "smooth_factor": 2,  # % 'Smoothing factor:'
    "im_smooth_cycles": 0,  # % 'Smoothing cycles image:'
    "phasor_smooth_cycles": 1,  # % 'Smoothing cycles phasor:'
    "foreground_threshold": 10,
    "tau_exc": np.inf,  # % 'Tau_exc'
    "intercept_at_origin": False,  # % 'Fix Intercept at origin'

    # Parameters that are calculated in th matlab code but that could be modified manually
    "M0": None,
    "Min": None,

    # Paramaters that are fixed in the matlab code
    "m0": 1,
    "harm1": 1,  # MATLAB code: harm1=1+2*(h1-1), where h1=1
    "klog": 4,
}
# -----------------------------------------------------------
#    Paths to folders containing the images, 1 per dye.
f1=easygui.diropenbox(default=os.path.expanduser("~Desktop"))
f2=easygui.diropenbox(default=os.path.expanduser("~Desktop"))

filenames = [f1,f2]

#    List of depletion powers to be analyzed, string to be identified in the file name
powers=["_30"]
powersnum=[0]

#    List of keys to be used to extract the images from the msr files
keys=['Confocal_561 {11}']
#keys=[0]


colors_Centroids=["k","k","k"]
colors=['gray', "hotpink"]
names= ['Bassoon CF594', 'Homer STAR Orange']
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

#For each dye, extract the images from the msr files and calculate the phasor distribution, their centroids and ellipses (70th percentile)
for k,filename in enumerate(filenames) :
    
        # For each depletion power, list the images acquired with this power 
    for a,power in enumerate(powers):
        # Make list of all the images acquired with a certain depletion power in the folder
        extension = ".msr"
        path = os.path.join(filename, '*{}PercentSTED.msr'.format(power) )
        images = glob.glob(path)
        print('There are ',len(images), ' msr files in this folder')
        if len(images) == 0:
            path = os.path.join(filename, '*{}PercentSTED.tiff'.format(power) )
            images = glob.glob(path)
            print('There are ',len(images), ' tiff files in this folder')
            extension = ".tiff"


        msrfiles=images

        #For each image, calculate the phasor distribution and plot the scatter plot of the phasor distribution
        for i, msr in enumerate(msrfiles) :
    # Read the msr file
            imagemsr = load_image(msr)
            print(os.path.basename(msr))
            #print(imagemsr.keys())

    # Calculate the phasor distribution for all pixels in the image that are above the foreground threshold
            df = pd.DataFrame(columns=['x','y'])
            dg = pd.DataFrame(columns=['g', 's'])
            image1=select_channel(imagemsr, keys[a])
            #image1=imagemsr[keys[a]]
            print("Caclulation for an image of shape", image1.shape, "...")
            #params_dict["foreground_threshold"] = get_foreground(image1)
            params_dict["foreground_threshold"]=10
            print("foreground_threshold=", params_dict["foreground_threshold"])
            x,y,g_smoothed,s_smoothed, orginal_idxs= Median_Phasor(image1, params_dict, **params_dict, show_plots=False)
            df['x']=x.flatten()
            df['y']=y.flatten()
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

            ax_centroids.scatter(dg['g'], dg['s'],s=0.5, c=colors[1],alpha=0.1,rasterized=True)
            calibx,caliby=polar_to_cart([IRF[0],1],[IRF[1],0])
            ax_centroids.scatter(calibx,caliby,color=colors_Centroids[1], s=25)
            

    # Calculate the centroid of the phasor distribution
            kmeans = KMeans(n_clusters=1, init='k-means++', random_state=42)
            y_kmeans = kmeans.fit_predict(dg)
            CoM_x=kmeans.cluster_centers_[:, 0][:]
            CoM_y=kmeans.cluster_centers_[:, 1][:]

            ax_centroids.scatter(CoM_x, CoM_y, color=colors_Centroids[2], s=25)


    


fig_centroids.savefig("Centroids_single.pdf",transparent=True, bbox_inches="tight",dpi=900)

plt.show()
