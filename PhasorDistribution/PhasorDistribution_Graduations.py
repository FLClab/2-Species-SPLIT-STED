
"""
Script to plot a phasor distribution before and after calibration with the IRF measurement
"""
import skimage

import os
import glob
import numpy 
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
import random
import scipy
from shapely.geometry.point import Point
from shapely import affinity
from sklearn.cluster import KMeans
matplotlib.rcParams['axes.linewidth'] = 0.8



colors_Centroids=["k","m","k"]
colors=['gray', "hotpink"]

# Phasor coordinates of the centroid of the IRF measurement in polar coordinates
IRF= (0.9527011687260826, 0.4695955819269703)
MeanPositions={}
Ellipsedims={}



#Create the figure and axes for the 2D plot of centroids and ellipses 
fig_centroids,ax_centroids = plt.subplots(figsize=(4,3))
theta = numpy.linspace(0, numpy.pi, 100)
r = 0.5
x1 = r * numpy.cos(theta) + 0.5
x2 = r * numpy.sin(theta)
ax_centroids.plot(x1, x2, color="black", ls="--",linewidth=0.8)
ax_centroids.set_xlim(-0.02, 1.02)
ax_centroids.set_ylim(-0.02, 0.82)
ax_centroids.set_xlabel('g')
ax_centroids.set_ylabel('s')
w=2*numpy.pi*40E6
taus=[0,1,2,3,4,5,6,7,8]
for tau in taus:
    g_tau = (1/(1+(w*tau*1E-9)**2))
    s_tau = (w*tau*1E-9/(1+(w*tau*1E-9)**2))
    ax_centroids.plot(g_tau, s_tau, 'o', markersize=8,label=str(tau)+' ns')
    #ax_centroids.text(g_tau+0.02, s_tau, str(tau)+' ns')

g_tau1=(1/(1+(w*5*1E-9)**2))
g_tau2=(1/(1+(w*1*1E-9)**2))
s_tau1=(w*5*1E-9/(1+(w*5*1E-9)**2))
s_tau2=(w*1*1E-9/(1+(w*1*1E-9)**2))
ax_centroids.plot(g_tau1/2+g_tau2/2, s_tau1/2+s_tau2/2, 'o', markersize=8,label='Mix 5ns & 1ns 1/2')
ax_centroids.plot(g_tau1/3+(2*g_tau2/3), s_tau1/3+(2*s_tau2/3), 'o', markersize=8,label='Mix 5ns & 1ns 1/3')
ax_centroids.plot(g_tau2/3+(2*g_tau1/3), s_tau2/3+(2*s_tau1/3), 'o', markersize=8,label='Mix 5ns & 1ns 1/3')
ax_centroids.legend()

fig_centroids.savefig("Phasor_Graduations.pdf",transparent=True, bbox_inches="tight")



values=numpy.linspace(0,250,num=250,endpoint=False)
tvalues=values*0.08
tvalues=tvalues+0.08
print(tvalues)
taus=[10,1,2,3,4,5,6,7,8]
fig,ax=plt.subplots(figsize=(4,3))
for tau in taus:
    distribution=numpy.exp(-tvalues/tau)
    #ax.plot(tvalues, distribution, label=f'tau={tau}',linewidth=2)
    tlist=random.choices(tvalues, weights=distribution, k=100000)
    #noisegrid=numpy.random.poisson(lam=l, size=10000).astype(int)
    ax.hist(tlist, bins=250, range=(0.08, 20),histtype='step',fill=False, label=f'tau={tau}',linewidth=2)
distribution_tau1=numpy.exp(-tvalues/5)
distribution_tau2=numpy.exp(-tvalues/1)
tlist_mix1=random.choices(tvalues, weights=distribution_tau1/2+distribution_tau2/2, k=100000)
tlist_mix2=random.choices(tvalues, weights=distribution_tau1/3+2*(distribution_tau2/3), k=100000)
tlist_mix3=random.choices(tvalues, weights=2*(distribution_tau1/3)+(distribution_tau2/3), k=100000)
ax.hist(tlist_mix1, bins=250, range=(0.08, 20),histtype='step', fill=False, label='Mix 5ns & 1ns',linewidth=2)
ax.hist(tlist_mix2, bins=250, range=(0.08, 20),histtype='step', fill=False, label='Mix 5ns & 1ns',linewidth=2)
ax.hist(tlist_mix3, bins=250, range=(0.08, 20),histtype='step', fill=False, label='Mix 5ns & 1ns',linewidth=2)

ax.legend()
ax.set_xlabel("Time (ns)")
fig.savefig("Histograms_graduations.pdf",transparent=True, bbox_inches="tight")
plt.show()
