
"""
This script calculates phasor distributions of 2 different dyes imaged independently
it then calculates the centroid of the phasor distribution and the ellipse that best fits the distribution

it then calculates for all pairs of images of the same power the following metrics: 
-the distance between the centroids 
-the intersection over union of the ellipses
-the shortest distance between the ellipses

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
f1=easygui.diropenbox(default=os.path.expanduser("~Desktop"), title="Select folder with images of first dye")
f2=easygui.diropenbox(default=os.path.expanduser("~Desktop"), title="Select folder with images of second dye")

# Ask user to input the name of the folder to save the results
savefoldername =str(input("Name of folder to save to: "))
  
filenames = [f1,f2]

# List of powers to use for the phasor calculation 
# string to look for in the filenames
powers=["_*","_5","_10","_15","_20"]
#powers=["_*","_10","_20","_30","_40"]

# Values to use for the powers in the plot
#powersnum=[0,10,20,30,40]
powersnum=[0,5,10,15,20]

ticklabels=["0","22","44","66","88"]
mwpowers=["","0","44","88"]

# List of channels to use for the phasor calculation. For Tiff files, use the channel number instead of the key

#keys=['Conf_ 594 {2}','STED_594 {2}','STED_594 {2}','STED_594 {2}','STED_594 {2}']
keys=['Conf_635P {2}', 'STED_635P {2}', 'STED_635P {2}', 'STED_635P {2}', 'STED_635P {2}', 'STED_635P {2}']
#keys=['Conf640 {10}','STED640 {10}','STED640 {10}','STED640 {10}','STED640 {10}']
keys=['Confocal_561 {11}','STED 561 {11}','STED 561 {11}','STED 561 {11}','STED 561 {11}']
#keys=['Conf_635P {2}','Conf_635P {2}','Conf_635P {2}','Conf_635P {2}']
#keys=['STED 561 {11}','STED 561 {11}','STED 561 {11}','STED 561 {11}']

keys=[0,1,1,1,1] # For Tiff files, use the channel number instead of the key

colors=['magenta','cyan','lawngreen','mediumpurple']
#colors_centroids=[['lightskyblue', 'deepskyblue','blue','mediumblue','darkblue',"cyan"],["pink","lightpink",'lightcoral','indianred','mediumvioletred',"magenta"]]
colors=[['deepskyblue', 'deepskyblue','deepskyblue','deepskyblue','deepskyblue','deepskyblue'],["hotpink","hotpink","hotpink","hotpink","hotpink","hotpink"]]
colors_centroids=[["#7ce8ffff","#55d0ffff","#00acdfff","#0080bfff","#00456bff"],["#fcbcd7ff","#f9a3cbff","#ef87beff","#e569b3ff","#bf4290ff"]]
#colors=["orangered"]

names= ['Bassoon CF594', 'Homer STAR Orange']


Positions={}
MeanPositions={}
Ellipsedims={}
Filenames={}
# Create a dataframe to store the results
Overall_data = pd.DataFrame(columns=["Power",'image1', 'image2', 'IOU','Centroid Distance','Shortest Distance'])

# Create a folder to save the results
savefolder=os.path.join(os.path.expanduser("~/Desktop"),savefoldername + "_PhasorDists")
os.makedirs(savefolder,exist_ok=True)
    

# Create a figure to plot the phasors
fig,ax_scatter= plt.subplots(figsize=(3,3))
ax_scatter.set_xlim(0, 1)
ax_scatter.set_ylim(0, 1)


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
            imagemsr = load_image(msr)
            print(os.path.basename(msr))
           
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
    

table=np.zeros((2,len(powersnum)))
table2=np.empty((2,1))
table4=np.empty((2,1))
table3=np.zeros((2,len(powersnum)))
table5=np.empty((2,1))
table6=np.zeros((2,len(powersnum)))
print(table2.shape)
start=0
# Loop through the powers and plot the centroids and ellipses for each dye
for b,power in enumerate(powersnum):
    print("power",power)
    # Plot the centroids of the phasor distributions
    ax_scatter.scatter(MeanPositions[0][power][0], MeanPositions[0][power][1], s=10, c=colors_centroids[0][b])
    ax_scatter.scatter(MeanPositions[1][power][0], MeanPositions[1][power][1], s=10, c=colors_centroids[1][b])
    # Create the ellipses for the phasor distribution of each dye and plot them
    ell1 = mpatches.Ellipse([MeanPositions[0][power][0], MeanPositions[0][power][1]], width=np.mean(Ellipsedims[0][power][0, :]),
                           height=np.mean(Ellipsedims[0][power][1, :]), angle=np.mean(Ellipsedims[0][power][2, :]))
    ell2 = mpatches.Ellipse([MeanPositions[1][power][0], MeanPositions[1][power][1]], width=np.mean(Ellipsedims[1][power][0, :]),
                           height=np.mean(Ellipsedims[1][power][1, :]), angle=np.mean(Ellipsedims[1][power][2, :]))
    
    ellipse1=create_ellipse((MeanPositions[0][power][0], MeanPositions[0][power][1]), (np.mean(Ellipsedims[0][power][0, :]),np.mean(Ellipsedims[0][power][1, :])), angle=np.mean(Ellipsedims[0][power][2, :]))
    ellipse2=create_ellipse((MeanPositions[1][power][0], MeanPositions[1][power][1]), (np.mean(Ellipsedims[1][power][0, :]),np.mean(Ellipsedims[1][power][1, :])), angle=np.mean(Ellipsedims[1][power][2, :]))
    intersect = ellipse1.intersection(ellipse2) # intersection of the two ellipses


    ax_scatter.add_artist(ell1)
    ax_scatter.add_artist(ell2)
    ell1.set_facecolor("None")
    ell2.set_facecolor("None")
    ell1.set_edgecolor(colors_centroids[0][b])
    ell2.set_edgecolor(colors_centroids[1][b])
    ell1.set_linewidth(0.8)
    ell2.set_linewidth(0.8)
    table[0,b]=power
    # Calculate the distance between the centroids of the phasor distributions
    table[1,b]=scipy.spatial.distance.euclidean([MeanPositions[0][power][0],MeanPositions[0][power][1]],[MeanPositions[1][power][0],MeanPositions[1][power][1]])
    table3[0,b]=power
    # Calculate the intersection over union (IoU) for the ellipses
    table3[1,b]=intersect.area/(ellipse1.area+ellipse2.area)
    table6[0, b] = power
    # Calculate the shortest distance between the ellipse of each dye
    table6[1, b] =ellipse1.distance(ellipse2)



    alldists=scipy.spatial.distance.cdist(Positions[0][power].T,Positions[1][power].T, metric='euclidean').ravel()
    #print(alldists)
    print(Positions[0][power].shape,Positions[1][power].shape)
    num = alldists.shape[0]
    print("num",num)
    tabletemp=np.zeros((2,num))
    tabletemp4 = np.zeros((2, num))
    tabletemp5 = np.zeros((2, num))
    for di,dist in enumerate(alldists):
        tabletemp[0, di] = power
        tabletemp[1, di] =dist
    
    # Calculate the intersection over union (IoU) and shortest distance between the ellipses for each possible pair of images (1 per dye)
    for d,pair in enumerate(itertools.product(range(Positions[0][power].shape[1]),range(Positions[1][power].shape[1]))):

        ellipse1 = create_ellipse((Positions[0][power][0,pair[0]], Positions[0][power][1,pair[0]]),
                                  (Ellipsedims[0][power][0, pair[0]], Ellipsedims[0][power][1, pair[0]]),
                                  angle=Ellipsedims[0][power][2, pair[0]])
        ellipse2 = create_ellipse((Positions[1][power][0,pair[1]], Positions[1][power][1,pair[1]]),
                                  (Ellipsedims[1][power][0, pair[1]], Ellipsedims[1][power][1, pair[1]]),
                                  angle=Ellipsedims[1][power][2, pair[1]])

        intersect = ellipse1.intersection(ellipse2)
     
        #IOU=intersect.area / (ellipse1.area + ellipse2.area-intersect.area)
        IOU=intersect.area / (ellipse1.union(ellipse2).area)
        tabletemp4[0, d] = power
        tabletemp5[0, d] = power
        tabletemp4[1, d] =IOU
        tabletemp5[1, d] = ellipse1.distance(ellipse2)
        #print(Overall_data)
        centroiddist=tabletemp[1,d]
        Overall_data.loc[Overall_data.shape[0]] = [power,Filenames[0][power][pair[0]],Filenames[1][power][pair[1]],IOU,centroiddist,ellipse1.distance(ellipse2)]
    table4 = np.concatenate([table4, tabletemp4], axis=1)
    table5 = np.concatenate([table5, tabletemp5], axis=1)


       


    table2=np.concatenate([table2,tabletemp],axis=1)
    #table2=np.concatenate([table2,tabletemp],axis=1)
    #print(table2.shape)
    start+=num

# Save the results in the dataframe to a csv file
Overall_data.to_csv(os.path.join(savefolder,"Overall_data_Ellipses_"+savefoldername+".csv"))
print(Overall_data)
table2=table2[:,1:]
table4=table4[:,1:]
table5=table5[:,1:]

# Calculate the mean and standard deviation of the distances between the centroids and the ellipses for each power
table4mean=np.empty((3,1))
listtable4=[]
for x in sorted(np.unique(table4[0,:])):
    a=table4[0,:]
    b = table4[1, :]
    temp=np.array([np.mean(b[np.where(a==x)]),np.std(b[np.where(a==x)]),x])
    table4mean = np.column_stack([table4mean, temp])
    listtable4.append(b[np.where(a==x)])

table5mean=np.empty((3,1))
listtable5=[]
for x in sorted(np.unique(table5[0,:])):
    a=table5[0,:]
    b = table5[1, :]
    temp=np.array([np.mean(b[np.where(a==x)]),np.std(b[np.where(a==x)]),x])
    listtable5.append(b[np.where(a==x)])

    table5mean = np.column_stack([table5mean, temp])
table4mean=table4mean[:,1:]
table5mean=table5mean[:,1:]
#print("table5mean",table5mean.shape,table5mean)
#print("table4mean",table4mean.shape,table4mean)
table2mean=np.empty((3,1))
listtable2=[]
for x in sorted(np.unique(table2[0,:])):
    a=table2[0,:]
    b = table2[1, :]
    temp=np.array([np.mean(b[np.where(a==x)]),np.std(b[np.where(a==x)]),x])
    listtable2.append(b[np.where(a==x)])
    table2mean = np.column_stack([table2mean, temp])

fig1,ax1=plt.subplots(figsize=(3,1.5))

ax1.scatter(table2[0,:],table2[1,:],c="blueviolet",s=25,rasterized=True)
#seaborn.stripplot(x=table2[0,:], y=table2[1,:],ax=ax1, size=2,color="blueviolet")
print("Unique",np.unique(table2[0,:]))
ax1.violinplot(listtable2,positions=powersnum,widths=8)
ax1.set_ylabel("Distance between centroids (different dye)")
ax1.set_xticks(ticks=powersnum,labels=ticklabels)
ax1.set_xlabel("Depletion Power [mW]")
ax1.set_title("All centroid distances")

fig11 = plt.figure(figsize=(12,8))
plt.scatter(table3[0,:],table3[1,:])
plt.title("Mean IOU")
fig2,ax2 = plt.subplots(figsize=(2,1))
fig111,ax111 = plt.subplots(figsize=(2,1))
#plt.bar(table6[0,:],table6[1,:],width=5)
ax111.violinplot(listtable5,positions=powersnum,widths=8)
ax2.violinplot(listtable5,positions=powersnum,widths=8)
#plt.bar(table5mean[2,:],table5mean[0,:],yerr=table5mean[1,:],width=5,lw=3,capsize=10)
ax111.scatter(table5[0,:],table5[1,:],c="blueviolet",s=15,rasterized=True)
#seaborn.stripplot(x=table5[0,:], y=table5[1,:],ax=ax111, size=2,color="blueviolet")
ax111.set_ylabel("Shortest distance between ellipses (different dye)")
ax2.set_ylabel("Shortest distance between ellipses (different dye)")
#plt.xticks(ticks=powersnum,labels=ticklabels)
ax111.set_xlabel("Depletion Power [mW]")
ax2.set_xlabel("Depletion Power [mW]")
fig1111 = plt.figure(figsize=(3,1.5))


## Compute the statistics to determine if the metrics are significantly different betweeen the different depletion powers

print("################################################################")
print("STATS IOU")
combine=list(itertools.chain.from_iterable(listtable4))
print(len(listtable4))
print(len(combine))
stat, p = scipy.stats.shapiro(combine)
print("Shapiro result",stat, p)
if p<0.05:
    print("All together,IOUs are not normal")

for t4 in listtable4:
    stat, p = scipy.stats.shapiro(t4)
    print("Shapiro result",stat, p)
    if p<0.05:
        print("IOUs are not normal")
try:
    stat, p =scipy.stats.kruskal(*listtable4)
    print("kruskal result",stat, p)
except ValueError:
    print("Kruskal did not work")
print("###############################################################")
combine=list(itertools.chain.from_iterable(listtable5))
print("STATS Shortest ditances")
stat, p = scipy.stats.shapiro(combine)
print("Shapiro result",stat, p)
if p<0.05:
    print("All together,Shortest distances are not normal")

for t4 in listtable4:
    stat, p = scipy.stats.shapiro(t4)
    print("Shapiro result",stat, p)
    if p<0.05:
        print("IOUs are not normal")
try:
    stat, p =scipy.stats.kruskal(*listtable4)
    print("kruskal result",stat, p)
except ValueError:
    print("Kruskal did not work")
    

sig=get_significance(listtable4, verbose=True)
print("sig",sig)


plt.violinplot(listtable4,positions=powersnum,widths=8)
plt.scatter(table4[0,:],table4[1,:],c="blueviolet",s=25,rasterized=True)
plt.ylabel("Intersection over union of ellipses")
plt.xticks(ticks=powersnum,labels=ticklabels)
plt.xlabel("Depletion Power [mW]")

# save the violin plot figures
fig111.savefig(os.path.join(savefolder,'Violinplot_Shortestdistance_Ellipses.pdf'),transparent='True', bbox_inches="tight",dpi=900)
fig1111.savefig(os.path.join(savefolder,'Violinplot_IOU_Ellipses.pdf'),transparent='True', bbox_inches="tight",dpi=900)
fig4,ax4 = plt.subplots(figsize=(6,10))
fig2,ax2 = plt.subplots(ncols=2,nrows=2,sharex=True,sharey='row',figsize=(15,12))
fig3,ax3 = plt.subplots(ncols=2,nrows=2,sharex=True,sharey='row',figsize=(15,12))
fig5,ax5=plt.subplots(nrows=2,ncols=1,figsize=(8,10),sharex=True)
ax2[1,0].set_xlabel("Depletion Power [mW]",fontsize=16)
ax2[1,1].set_xlabel("Depletion Power [mW]",fontsize=16)
ax3[1,0].set_xlabel("Depletion Power [mW]",fontsize=16)
ax3[1,1].set_xlabel("Depletion Power [mW]",fontsize=16)
ax4.set_xlabel("Depletion Power [mW]",fontsize=16)
ax5[1].set_xlabel("Depletion Power [mW]",fontsize=16)
ax5[1].set_xlim(-0.1, 1)
ax5[1].set_ylim(0.4, 0.6)
ax5[0].set_xlim(-0.1, 1)
ax5[0].set_ylim(0.4, 0.6)
ax2[0,0].set_ylabel("Ellipse Minor Axis length",fontsize=16)
ax2[1,0].set_ylabel("Ellipse Major Axis Length ",fontsize=16)
ax3[0,0].set_ylabel("Ellipse Aspect Ratio",fontsize=16)
ax3[1,0].set_ylabel("Ellipse Area ",fontsize=16)

positions=[0,0.2,0.4,0.6,0.8,1.0]

ax5[1].set_xticklabels(mwpowers,fontsize=16)
posdif=[-1,1]
colorss=["deepskyblue","hotpink"]

lines = [mpatches.Patch(color=colorss[i], label=names[i]) for i in range(len(names))]

ax2[0,0].legend(handles=lines)


ax2[1,0].set_xticklabels(mwpowers,fontsize=16)
ax2[1,1].set_xticklabels(mwpowers,fontsize=16)
ax3[1,0].set_xticklabels(mwpowers,fontsize=16)
ax3[1,1].set_xticklabels(mwpowers,fontsize=16)

for k,key in enumerate(Ellipsedims.keys()):

    for p,power in enumerate(Ellipsedims[key].keys()):

        x=power*np.ones((1, Ellipsedims[key][power].shape[1]))
        semimajor=(Ellipsedims[key][power][1,:]/2)
        semiminor=(Ellipsedims[key][power][0,:]/2)

        ecc = semiminor / semimajor
        area=np.pi*semiminor*semimajor

        ell = mpatches.Ellipse([positions[p],0.5],width=np.mean(Ellipsedims[key][power][0,:]),height=np.mean(Ellipsedims[key][power][1,:]),angle=np.mean(Ellipsedims[key][power][2,:]))

        ax5[k].add_artist(ell)
        ell.set_facecolor("None")
        ell.set_edgecolor(colors[k][p])
        ell.set_linewidth(3)


        ax2[0,k].bar(power,np.mean(Ellipsedims[key][power][0,:]),yerr=np.std(Ellipsedims[key][power][0,:]),width=5,lw=3,capsize=10,color=colors[k][p])
        ax2[1,k].bar(power,np.mean(Ellipsedims[key][power][1,:]),yerr=np.std(Ellipsedims[key][power][1,:]),width=5,lw=3,capsize=10,color=colors[k][p])
        ax3[0,k].bar(power,np.mean(ecc),yerr=np.std(ecc),color=colors[k][p],width=5,lw=3,capsize=10)
        ax3[1, k].bar(power, np.mean(area), yerr=np.std(area), color=colors[k][p], width=5, lw=3, capsize=10)
        ax2[0,k].scatter(x,Ellipsedims[key][power][0,:],c=colors[k][p],edgecolors="black")
        ax2[1,k].scatter(x,Ellipsedims[key][power][1,:],c=colors[k][p],edgecolors="black")
        ax3[0,k].scatter(x,ecc,c=colors[k][p],edgecolors="black")
        ax3[1, k].scatter(x, area, c=colors[k][p], edgecolors="black")

    
fig2.savefig(os.path.join(savefolder,'Barplot_MajorMinorlengths_Ellipses.pdf'),transparent='True', bbox_inches="tight")
fig3.savefig(os.path.join(savefolder,'Barplot_AreaAspectratio_Ellipses.pdf'),transparent='True', bbox_inches="tight")
fig5.savefig(os.path.join(savefolder,'Mean_Ellipses.pdf'),transparent='True', bbox_inches="tight")
fig.savefig(os.path.join(savefolder,"Centroids_Mean.pdf"),transparent=True, bbox_inches="tight",dpi=900)
fig1.savefig(os.path.join(savefolder,"Centroids_Distances_violin.pdf"),transparent=True, bbox_inches="tight",dpi=900)
plt.close('all')
