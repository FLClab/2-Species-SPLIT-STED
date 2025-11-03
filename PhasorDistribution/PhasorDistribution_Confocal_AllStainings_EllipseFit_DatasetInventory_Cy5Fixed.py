
"""
This script calculates phasor distributions of a list of different dyes imaged independently.
It then calculates the centroid of the phasor distribution and the ellipse that best fits the 70th percentile of the distribution

It then plots the centroid and ellipse for each dye on the same phasor plot

and saves the results in a csv file
"""
import skimage
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn
import os
import itertools
import easygui
import math
from sys import path as path1; 
Functionspath=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Functions")
path1.append(Functionspath)
from statistics_functions import get_significance
from Main_functions import (load_image,select_channel,to_polar_coord, polar_to_cart, get_foreground)
from Phasor_functions import Median_Phasor
import scipy
from shapely.geometry.point import Point
from shapely import affinity
from sklearn.cluster import KMeans
matplotlib.rcParams['axes.linewidth'] = 0.8
#plt.style.use('dark_background')

# Functions to create ellipses 
def create_ellipse(center, lengths, angle=0):
    """
    create a shapely ellipse. adapted from
    https://gis.stackexchange.com/a/243462
    """
    circ = Point(center).buffer(1)
    #print(circ.area)
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

# ------------------ Default Input variables ----------------
params_dict = {
    
    "smooth_factor": 0.2,  # % 'Smoothing factor:'
    "phasor_smooth_cycles": 1,  # % 'Smoothing cycles phasor:'
    "foreground_threshold": 10,
    "harm1": 1,
}

# Path to folder containing the images, 1 per dye. A file browser will open for the user to navigate to and select the desired folder.

# f1 = os.path.join('U:', os.sep,'adeschenes','2023-12-04_FLIM_SynapticProteins_MediumAcquisition','Bassoon_CF594_STEDPowerBleach_MediumAcq_1')
# f2 = os.path.join('U:', os.sep,'adeschenes','2023-12-04_FLIM_SynapticProteins_MediumAcquisition','Homer_STOrange_STEDPowerBleach_MediumAcq_1')
# f3 = os.path.join('U:', os.sep,'adeschenes',"2023-12-21_FLIM_MediumAcq_Spectrin_Actin_Bassoon","B2Spectrin_CF594_STEDPowerBleach_MediumAcq_1")
# f4 = os.path.join('U:', os.sep,'adeschenes',"2023-12-21_FLIM_MediumAcq_Spectrin_Actin_Bassoon","B2Spectrin_STOrange_STEDPowerBleach_MediumAcq_1")
# f5 = os.path.join('U:', os.sep,'adeschenes',"2023-12-21_FLIM_MediumAcq_Spectrin_Actin_Bassoon","Bassoon_CF594_STEDPowerBleach_MediumAcq_1")
# f6=os.path.join('U:', os.sep,'adeschenes',"2024-03-06_FLIM_PSDBassoon_Cy3","msPSD95_STOrange_STEDPowerBleach_MediumAcq_MoreReps_1")
# f7=os.path.join('U:', os.sep,'adeschenes',"2024-03-06_FLIM_PSDBassoon_Cy3","rabBassoon_CF594_STEDPowerBleach_MediumAcq_MoreReps_1")




f1 = os.path.join('U:', os.sep,'mlafontaine',"22-07-06_MIXEDFLIM","Bassoon_ATTO647N-2")

f2 = os.path.join('U:', os.sep,'mlafontaine',"22-07-06_MIXEDFLIM","Tubuline_STAR635p-5")


###############################################################################################################################



f5 = os.path.join('U:', os.sep,'adeschenes',"2022-12-09_STEDFLIM_BassoonHomer","Cy5","STEDPower_Bleach_DifferentRegions","Bassoon_STAR635P_GAM_STEDPower_Bleach_1")
f6 = os.path.join('U:', os.sep,'adeschenes',"2022-12-09_STEDFLIM_BassoonHomer","Cy5","STEDPower_Bleach_DifferentRegions","Bassoon_STAR635P_Nanotag_STEDPower_Bleach_1")
f7 = os.path.join('U:', os.sep,'adeschenes',"2022-12-09_STEDFLIM_BassoonHomer","Cy5","STEDPower_Bleach_DifferentRegions","Homer_Atto647N_STEDPower_Bleach_1")


f8= os.path.join('U:', os.sep,'adeschenes',"2023-02-16_FLIM_SynapticProteins","Bassoon_Atto647N_STEDPower_Bleach_60PControl_1")
f9= os.path.join('U:', os.sep,'adeschenes',"2023-02-16_FLIM_SynapticProteins","PSD95_ST635P_STEDPower_Bleach_60PControl_1")
f10= os.path.join('U:', os.sep,'adeschenes',"2023-02-16_FLIM_SynapticProteins","PSD95_STRed_STEDPower_Bleach_60PControl_2")


f11= os.path.join('U:', os.sep,'adeschenes',"2023-02-19_FLIM_SynapticProteins_2","Bassoon_Atto647N_STEDPower_5to40_Bleach_60P_Control_1")
f12= os.path.join('U:', os.sep,'adeschenes',"2023-02-19_FLIM_SynapticProteins_2","PSD95_ST635P_STEDPower_5to40_Bleach_60P_Control_1")
f13= os.path.join('U:', os.sep,'adeschenes',"2023-02-19_FLIM_SynapticProteins_2","PSD95_STRED_STEDPower_5to40_Bleach_60P_Control_1")

                  
f14= os.path.join('U:', os.sep,'adeschenes',"2023-03-10_FLIM_SynaticProteins_3","Bassoon_Atto647N_STEDPowerBleach_5to40_60PControl_1")

f15= os.path.join('U:', os.sep,'adeschenes',"2023-03-10_FLIM_SynaticProteins_3","PSD95_ST635P_STEDPowerBleach_5to40_60PControl_1")


f16= os.path.join('U:', os.sep,'adeschenes',"2023-04-25_CalgarySampleTests","FLIM","Spectrin_Atto647N_STEDPowerBleach_5to40_60PControl")
f17= os.path.join('U:', os.sep,'adeschenes',"2023-04-25_CalgarySampleTests","FLIM","Phalloidin_ST635_STEDPowerBleach_5to40_60PControl")


f18= os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5',"PSD95_AF647_STEDPowerBleach_5to30_1")
f19=os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5','rab_Bassoon_STAR635P_STEDPowerBleach_5to30_1')
f20= os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5',"msB2Spectrin_AF647_STEDPowerBleach_5to30_1")
f21= os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5',"alphaTubulin_AF647_STEDPowerBleach_5to20_1")

f22=os.path.join('U:', os.sep,'adeschenes',"2024-03-06_FLIM_PSDBassoon_Cy3","AlphaTubulin_AF647_1to250_STEDPowerBleach_1")
f23=os.path.join('U:', os.sep,'adeschenes',"2024-03-06_FLIM_PSDBassoon_Cy3","AlphaTubulin_AF647_1to500_STEDPowerBleach_1")
filenames = [f1,f2,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,f20,f21,f22,f23]





#f1= os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5',"PSD95_AF647_STEDPowerBleach_5to30_1")
#f2=os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5','rab_Bassoon_STAR635P_STEDPowerBleach_5to30_1')
#f3= os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5',"msB2Spectrin_AF647_STEDPowerBleach_5to30_1")
#f4= os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5',"alphaTubulin_AF647_STEDPowerBleach_5to20_1")

#f1= os.path.join('U:', os.sep,'adeschenes',"2023-04-25_CalgarySampleTests","FLIM","Spectrin_Atto647N_STEDPowerBleach_5to40_60PControl")
#f2= os.path.join('U:', os.sep,'adeschenes',"2023-04-25_CalgarySampleTests","FLIM","Phalloidin_ST635_STEDPowerBleach_5to40_60PControl")


# Ask user to input the name of the folder to save the results
savefoldername =str(input("Name of folder to save to: "))
  
#filenames = [f1,f2,f3,f4,f5,f6,f7]
#filenames = [f1,f2,f3,f4]
# Names of the dyes imaged in each folder, in the same order as the filenames list
names= ['PSD95 AF647', 'rabBassoon STAR635P','B2Spectrin AF647','alphaTubulin AF647']
names=[os.path.basename(f) for f in filenames]
names=[name.split('_STEDPower')[0] for name in names]
#names= ['Bassoon CF594', 'Homer STAR Orange','B2Spectrin CF594','B2Spectrin STOrange','Bassoon CF594','msPSD95 STOrange','rabBassoon CF594']
# List of powers to use for the phasor calculation
# string to look for in the filenames
powers=["_*"]
#powers=["_*","_10","_20","_30","_40"]

# Values to use for the powers in the plot
#powersnum=[0,10,20,30,40]
powersnum=[0]

ticklabels=["0","22","44","66","88"]
mwpowers=["","0","44","88"]

# List of channels to use for the phasor calculation. For Tiff files, use the channel number instead of the key


#keys=['Conf640 {10}']
#keys=['Confocal_561 {11}']
keys=['Conf_635P {2}']


#keys=[0,1,1,1,1] # For Tiff files, use the channel number instead of the key

colors=['magenta','cyan','lawngreen','mediumpurple']
#colors_centroids=[['lightskyblue', 'deepskyblue','blue','mediumblue','darkblue',"cyan"],["pink","lightpink",'lightcoral','indianred','mediumvioletred',"magenta"]]
colors=[['deepskyblue', 'deepskyblue','deepskyblue','deepskyblue','deepskyblue','deepskyblue'],["hotpink","hotpink","hotpink","hotpink","hotpink","hotpink"]]
colors_centroids=[["#55d0ffff"],["#0080bfff"],["#00456bff"],["#ef87beff"],["#e569b3ff"],["#bf4290ff"],["#800055ff"]]
#colors=["orangered"]




Positions={}
MeanPositions={}
Ellipsedims={}
Filenames={}


# Create a folder to save the results
savefolder=os.path.join(os.path.expanduser("~/Desktop"),savefoldername + "_PhasorDists")
os.makedirs(savefolder,exist_ok=True)
    

# Create a figure to plot the phasors
fig,ax_scatter= plt.subplots(figsize=(3,3))
ax_scatter.set_xlim(0, 1)
ax_scatter.set_ylim(0, 1)
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.tab20b.colors)

# Create the universal semi-circle and plot it
theta = np.linspace(0, np.pi, 100)
r = 0.5
x1 = r * np.cos(theta) + 0.5
x2 = r * np.sin(theta)
ax_scatter.plot(x1, x2, color="black", ls="--",linewidth=0.8)
ax_scatter.set_xlabel('g')
ax_scatter.set_ylabel('s')


#axcentroids.plot(x1, x2, color="black", ls="--")
for k,filename in enumerate(filenames) :
    Filenames[k]={}
    Positions[k]={}
    MeanPositions[k]={}
    Ellipsedims[k]={}
    
    for a,power in enumerate(powers):
        ## For each power in the list of powers find all the images in the folder that have the power in their name
        
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
        Positions[k][powersnum[a]]=np.zeros((2,len(msrfiles)))
        Ellipsedims[k][powersnum[a]]=np.zeros((3,len(msrfiles)))
        Filenames[k][powersnum[a]]=[]
        for i, msr in enumerate(msrfiles) :
        # Read the image and calculate its phasor distribution with median filtering

            print(os.path.basename(msr))
            try:
                imagemsr = load_image(msr)
            except RuntimeError:
                print("Could not load image", msr)
                continue
           
            Filenames[k][powersnum[a]].append(os.path.basename(msr))


            df = pd.DataFrame(columns=['x','y'])
            dg = pd.DataFrame(columns=['g', 's'])
            image1=select_channel(imagemsr, keys[a])
      
            print("Caclulation for an image of shape", image1.shape, "...")
            #params_dict["foreground_threshold"] = get_foreground(image1)
            params_dict["foreground_threshold"]=10
            print("foreground_threshold=", params_dict["foreground_threshold"])
            x,y,g_smoothed,s_smoothed, orginal_idxs= Median_Phasor(image1, params_dict, **params_dict)
            df['x']=x.flatten()
            df['y']=y.flatten()
        # Calibrate the phasor distribution in polar coordinates based on the IRF measurement and return to (g,s) coordinates
            m, phi = to_polar_coord(df['x'], df['y'])
            g,s =polar_to_cart(m, phi)
            dg['g'], dg['s'] = g, s

        # Calculate the centroid of the phasor distribution using KMeans clustering
            kmeans = KMeans(n_clusters=1, init='k-means++', random_state=42)
            y_kmeans = kmeans.fit_predict(dg)
            
            CoM_x=kmeans.cluster_centers_[:, 0][:]
            CoM_y=kmeans.cluster_centers_[:, 1][:]
           
            #print(a)
        # Calculate the ellipse that best fits the phasor distribution's 70th percentile
            cov = np.cov(g, s)
            val, vec = np.linalg.eig(cov)
            order = val.argsort()[::-1]

            eigen_val=val[order]
            norm_eigen_vec = vec[:,order]
            eigen_val = np.sort(np.sqrt(eigen_val))
            ppfs=[0.3]
            for ppf in ppfs:
                width = 2 * eigen_val[0] * np.sqrt(scipy.stats.chi2.ppf(ppf, 2))
                height = 2 * eigen_val[1] * np.sqrt(scipy.stats.chi2.ppf(ppf, 2))
                angle = np.rad2deg(np.arctan2(norm_eigen_vec[1, eigen_val.argmax()],
                                            norm_eigen_vec[0, eigen_val.argmax()]))
                ell = mpatches.Ellipse(dg.mean(axis=0),width=width,height=height,angle=angle)

            Positions[k][powersnum[a]][0,i]=CoM_x[0]

            Positions[k][powersnum[a]][1, i] =CoM_y[0]

            Ellipsedims[k][powersnum[a]][0, i] =width
            Ellipsedims[k][powersnum[a]][1, i] =height
            if angle >0:
                Ellipsedims[k][powersnum[a]][2, i] =angle
            else:
                Ellipsedims[k][powersnum[a]][2, i] = 180+angle

        MeanPositions[k][powersnum[a]]=np.mean(Positions[k][powersnum[a]],axis=1)
    

start=0
power=0
# Loop through the powers and plot the centroids and ellipses for each dye
for k,folder in enumerate(filenames):
    print("folder",power)
    # Plot the centroids of the phasor distributions
    ax_scatter.scatter(MeanPositions[k][power][0], MeanPositions[k][power][1], s=10, label=names[k])
    
    # Create the ellipses for the phasor distribution of each dye and plot them
    ell1 = mpatches.Ellipse([MeanPositions[k][power][0], MeanPositions[k][power][1]], width=np.mean(Ellipsedims[k][power][0, :]),
                           height=np.mean(Ellipsedims[k][power][1, :]), angle=np.mean(Ellipsedims[k][power][2, :]))

    


    ax_scatter.add_artist(ell1)

    ell1.set_facecolor("None")

    ell1.set_edgecolor("k")

    ell1.set_linewidth(0.8)

ax_scatter.legend(loc='upper right')  
fig.savefig(os.path.join(savefolder,"Centroids_FixedCy5.pdf"),transparent=True, bbox_inches="tight",dpi=900)

plt.close('all')
