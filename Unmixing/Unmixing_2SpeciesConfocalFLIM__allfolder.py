# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file. 
"""
from cgitb import handler
from re import T
#from termios import FF1
from turtle import color
from xmlrpc.client import Marshaller

import os
import glob
import numpy
import easygui
import skimage.io as skio
import numpy as np
from scipy import signal
from scipy.optimize import curve_fit
from scipy.stats import expon,exponnorm
from pandas.plotting import table
import pandas as pd
from sys import path as syspath;
#syspath.append('/Users/marielafontaine/Documents/GitHub/TCSPC/Analyse/')
syspath.append(os.path.join(os.path.dirname(os.path.dirname(syspath[0])),"TCSPC","Analyse"))
print(os.path.join(os.path.dirname(os.path.dirname(syspath[0])),"TCSPC","Analyse"))

import scipy
import time
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.gridspec import GridSpec
import skimage
import os
import easygui
import plotly.express as px #mettre matplotly
import tifffile
from plotly.offline import plot
import matplotlib.patches as mpatches
t0 = time.time()
from sklearn.cluster import KMeans

from tiffwrapper import LifetimeOverlayer,imsave
import os.path
from sys import path as path1;
#path1.append('/Users/marielafontaine/Documents/GitHub/Abberior-STED-FLIM/Functions')
dossier = os.path.expanduser("~/Documents/Github/Abberior-STED-FLIM/Functions")
path1.append(dossier)
from Main_functions import (line_equation, to_polar_coord, polar_to_cart, choose_msr_file, get_foreground)
from Phasor_functions import Median_Phasor,unmix2species
#from lifetime import *
from convertmsr_bioformatsAB import MSRReader
#plt.style.use('dark_background')
# ------------------ Default Input variables ----------------
params_dict = {
    # Parameter in option in the matlab code
    #    "Tg" : 6, #% 'First frame to sum:'
    "Nb_to_sum": 250,  # The Tg infered from this variable override Tg
    "smooth_factor":0.2,  # % 'Smoothing factor:'
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

matplotlib.rcParams['axes.linewidth'] = 0.8
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    name="coolspring",
    colors=["#00ffffff", "#ff00ffff"]
)
matplotlib.colormaps.register(cmap=cmap, force=True)
matplotlib.colormaps.register(cmap=cmap.reversed(), force=True)

# -----------------------------------------------------------

#    Sélection des images dans un même fichier avec easygui

#f1=easygui.diropenbox(default=os.path.expanduser("~Desktop"))
#f2=easygui.diropenbox(default=os.path.expanduser("~Desktop"))
#f3=easygui.diropenbox(default=os.path.expanduser("~Desktop"))

 


#f1=os.path.join('T:', os.sep,'adeschenes',"SimulationDataset_STEDFLIM","Cy3","Bassoon_CF594","MediumAcq")
#f2=os.path.join('T:', os.sep,'adeschenes',"SimulationDataset_STEDFLIM","Cy3","Homer_STORANGE","MediumAcq")
#f3=os.path.join('T:', os.sep,'adeschenes',"Dataset_Mixed_Images_Cy3","Homer_STOrange_Bassoon_CF594","MediumAcq")

#f1= os.path.join('U:', os.sep,'adeschenes',"2023-12-21_FLIM_MediumAcq_Spectrin_Actin_Bassoon","B2Spectrin_CF594_STEDPowerBleach_MediumAcq_1")
#f2= os.path.join('U:', os.sep,'adeschenes',"2023-12-21_FLIM_MediumAcq_Spectrin_Actin_Bassoon","Phalloidin_AF594_STEDPowerBleach_MediumAcq_1")
#f3= os.path.join('U:', os.sep,'adeschenes',"2023-12-21_FLIM_MediumAcq_Spectrin_Actin_Bassoon","Phalloidin_AF594_B2Spectrin_CF594_STEDPowerBleach_MediumAcq_3")
f1= os.path.join('U:', os.sep,'adeschenes','2024-03-06_FLIM_PSDBassoon_Cy3',"rabBassoon_CF594_STEDPowerBleach_MediumAcq_MoreReps_1")
f2= os.path.join('U:', os.sep,'adeschenes','2024-03-06_FLIM_PSDBassoon_Cy3',"msPSD95_STOrange_STEDPowerBleach_MediumAcq_MoreReps_1")
#f3= os.path.join('U:', os.sep,'adeschenes','2024-03-06_FLIM_PSDBassoon_Cy3',"msPSD95_STOrange_rabBassoon_CF594_STEDPowerBleach_MediumAcq_MoreReps_1")
f3=os.path.join('T:', os.sep,'adeschenes',"Dataset_Mixed_Images_Cy3","PSD95_STOrange_rabBassoon_CF594")
f2= os.path.join('U:', os.sep,'adeschenes',"2023-12-21_FLIM_MediumAcq_Spectrin_Actin_Bassoon","Bassoon_CF594_STEDPowerBleach_MediumAcq_1")
f1= os.path.join('U:', os.sep,'adeschenes',"2023-12-21_FLIM_MediumAcq_Spectrin_Actin_Bassoon","Phalloidin_AF594_STEDPowerBleach_MediumAcq_1")
f3= os.path.join('U:', os.sep,'adeschenes',"2023-12-21_FLIM_MediumAcq_Spectrin_Actin_Bassoon","Phalloidin_AF594_Bassoon_CF594_STEDPowerBleach_MediumAcq_1")
f2=os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5','rab_Bassoon_STAR635P_STEDPowerBleach_5to30_1')
f1= os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5',"msB2Spectrin_AF647_STEDPowerBleach_5to30_1")
f3=os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5',"msB2Spectrin_AF647_rabBassoon_STAR635P_STEDPowerBleach_5to30_1")
f1= os.path.join('U:', os.sep,'adeschenes',"2023-12-21_FLIM_MediumAcq_Spectrin_Actin_Bassoon","Bassoon_CF594_STEDPowerBleach_MediumAcq_1")
f2= os.path.join('U:', os.sep,'adeschenes',"2023-12-21_FLIM_MediumAcq_Spectrin_Actin_Bassoon","B2Spectrin_STOrange_STEDPowerBleach_MediumAcq_1")
f3= os.path.join('U:', os.sep,'adeschenes',"2023-12-21_FLIM_MediumAcq_Spectrin_Actin_Bassoon","Spectrin_STOrange_Bassoon_CF594_STEDPowerBleach_MediumAcq_1")
#f1=easygui.diropenbox(default=os.path.expanduser("~Desktop"))
#f2=easygui.diropenbox(default=os.path.expanduser("~Desktop"))
#f3=easygui.diropenbox(default=os.path.expanduser("~Desktop"))

f3=os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5',"alphaTubulin_AF647_Bassoon_STAR635P_STEDPowerBleach_5to20_1")
f2=os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5','rab_Bassoon_STAR635P_STEDPowerBleach_5to30_1')
f1=os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5',"alphaTubulin_AF647_STEDPowerBleach_5to20_1")

mapcomp = { 'Conf FLIM' : 'Confocal_561 {11}',
            'STED FLIM' :'STED 561 {11}' }


colors=['magenta',  'c' ,'lightgreen']
#labels = ['Bassoon STARORANGE', 'MAP2 CF594', 'Mélange']
labels = ['Bassoon CF594', 'PSD95 STAR Orange', 'Mélange']
filenamescontrol = [f1, f2]
filenamemixed=f3
keys = ['Confocal_561 {11}', 'Confocal_561 {11}','Confocal_561 {11}']
keys = [ 'Conf_635P {2}','Conf_635P {2}','Conf_635P {2}' ]
msrfiles = []

# -----------------------------------------------------------
#     Choisir un fichier msr parmi plusieurs
savefolder=str(input("Name of Output folder: "))
savefolder = os.path.join(os.path.expanduser("~/Desktop"), "Unmixing_"+savefolder+"_Confocal")
os.makedirs(savefolder, exist_ok=True)
with MSRReader() as msrreader:
    for filename in filenamescontrol :
        path = os.path.join(filename, '*.msr' )
        images = glob.glob(path)
        print('There are ',len(images), 'Images in this folder')
        for imagei in images:
            print(os.path.basename(imagei)) 
        numim = int(input('Fichier msr a extraire (1er=0): '))
        image = images[numim]
        msrfiles.append(image)
    print(msrfiles)

    path = os.path.join(filenamemixed, '*.msr')
    mixedimages = glob.glob(path)

#plt.style.use('dark_background')
    fig4,ax_scatter = plt.subplots(figsize=(2,2))


    ax_scatter.set_xlim(0, 1.05)
    ax_scatter.set_ylim(-0.05, 1)

    ax_scatter.set_xlabel('g')
    ax_scatter.set_ylabel('s')
    edge = np.linspace(start=0, stop=15, num=200)
    theta = np.linspace(0, np.pi, 100)
    r = 0.5
    x1 = r*np.cos(theta) + 0.5
    x2 = r*np.sin(theta)
    ax_scatter.plot(x1,x2, color = "black", ls = "--",linewidth=0.8)
    CoM_x, CoM_y = [], []
    d_melange = pd.DataFrame(columns=['g', 's'])
    with open(os.path.join(savefolder,'legend.txt'),'w') as data:
        data.write("Controls\n")

    scatterlist = []
    for i, msr in enumerate(msrfiles) : 
        df = pd.DataFrame(columns=['x','y'])
        dg = pd.DataFrame(columns=['g', 's'])
        with open(os.path.join(savefolder,'legend.txt'),'a') as data:
            data.write("{}\t{}\t{}\n".format(labels[i],keys[i],msr))
        imagemsr=msrreader.read(msr)
        image1 = imagemsr[keys[i]]
        print(image1.shape)
        #image1 =image1[10: -10, 10: -10,:]
        print(image1.shape)
        imsum = image1.sum(axis=2)
        imsum = imsum.astype('int16')
        
        seuil = get_foreground(imsum)

        print("Caclulation for an image of shape", image1.shape, "...")

        params_dict["foreground_threshold"] = seuil

        params_dict["Nb_to_sum"] = image1.shape[2]
        print("foreground_threshold=", params_dict["foreground_threshold"])
        
        x,y,g_smoothed,s_smoothed, original_idxes= Median_Phasor(image1, params_dict, **params_dict, show_plots=False)
        df['x']=x.flatten()
        df['y']=y.flatten()
        m, phi = to_polar_coord(df['x'], df['y'])
        g,s =polar_to_cart(m, phi)
        dg['g'], dg['s'] = g, s

        kmeans = KMeans(n_clusters = 1, init = 'k-means++', random_state = 42)
        y_kmeans = kmeans.fit_predict(dg)
        CoM_x.extend(kmeans.cluster_centers_[:, 0][:].tolist())
        CoM_y.extend(kmeans.cluster_centers_[:, 1][:].tolist())

        #ax_scatter.scatter(g,s, c=colors[i], alpha=0.01)
        a=ax_scatter.scatter(g, s, s=1, c=colors[i], alpha=0.1,label=labels[i],rasterized=True)
        scatterlist.append(a)
        #ax_scatter.scatter(kmeans.cluster_centers_[:, 0], kmeans.cluster_centers_[:, 1], c='black')


    xaxis = np.linspace(0, 1, 100)

    ##Calcul de Pn
    # Pn_x, Pn_y = CoM_x[0], CoM_y[0]
    norm = numpy.sqrt((CoM_x[0] - 0.5) ** 2 + (CoM_y[0] ** 2))
    Pn_x = 0.5 + (r * (CoM_x[0] - 0.5) / norm)
    Pn_y = 0 + r * (CoM_y[0] - 0) / norm
    P_n = numpy.array([Pn_x, Pn_y])
    norm = numpy.sqrt((CoM_x[1] - 0.5) ** 2 + (CoM_y[1] ** 2))
    P2_x = 0.5 + (r * (CoM_x[1] - 0.5) / norm)
    P2_y = 0 + r * (CoM_y[1] - 0) / norm
    p2 = numpy.array([P2_x, P2_y])

    ## Droite entre Pn - p2
    x2, y2 = [Pn_x, P2_x], [Pn_y,P2_y]
    m2, c2 = line_equation(x2, y2)
    y2 = m2 * xaxis + c2
    ax_scatter.plot(xaxis, y2, 'dodgerblue')
    pnscatter=ax_scatter.scatter(Pn_x,Pn_y, c='darkred' )
    p2scatter=ax_scatter.scatter(P2_x,P2_y, c='darkred' )
    fig4.savefig(os.path.join(savefolder,
                              "Phasor_2species_Confocal_ControlsOnly.pdf"),
                 transparent='True',
                 bbox_inches="tight", dpi=900)
    pnscatter.remove()
    p2scatter.remove()
    ax_scatter.lines[-1].remove()
    t = [scatter.remove() for scatter in scatterlist]
    for m,mixedimage in enumerate(mixedimages):
        print("***********************************************************")
        print("Working on image number ",m," out of ",len(mixedimages))
        print("***********************************************************")
        d_melange = pd.DataFrame(columns=['g', 's'])
        df = pd.DataFrame(columns=['x', 'y'])
        dg = pd.DataFrame(columns=['g', 's'])
        imagemsr = msrreader.read(mixedimage)
        image1 = imagemsr[keys[2]]
        print(image1.shape)

        imsum = image1.sum(axis=2)
        imsum = imsum.astype('int16')
     # Pour l'image du mélange on ne filtre pas le foreground pour bien classifier les pixels moins brillants
        params_dict["foreground_threshold"] = 5
        
        params_dict["Nb_to_sum"] = image1.shape[2]
        print("foreground_threshold=", params_dict["foreground_threshold"])

        x, y, g_smoothed, s_smoothed, original_idxes = Median_Phasor(image1, params_dict, **params_dict, show_plots=False)
        df['x'] = x.flatten()
        df['y'] = y.flatten()
        m, phi = to_polar_coord(df['x'], df['y'])
        g, s = polar_to_cart(m, phi)
        dg['g'], dg['s'] = g, s
        #ax.imshow(imsum, cmap='hot')
        d_melange['g'], d_melange['s'] = g, s

        kmeans = KMeans(n_clusters=1, init='k-means++', random_state=42)
        y_kmeans = kmeans.fit_predict(dg)
        CoM_x.extend(kmeans.cluster_centers_[:, 0][:].tolist())
        CoM_y.extend(kmeans.cluster_centers_[:, 1][:].tolist())

        p3 = d_melange[['g', 's']].to_numpy() #phaseur qui sera projeté
        imsum_flat_lin1, imsum_flat_lin2, Solve = unmix2species(p3, original_idxes, image1, P_n, p2)
        fraction1=imsum_flat_lin2.copy()
        fraction2 = imsum_flat_lin1.copy()
        imsum_flat_lin1*= imsum
        imsum_flat_lin2*= imsum

        mixphasor = ax_scatter.scatter(g, s, s=1, c=Solve[1,:],cmap="cool", label="Mixture",rasterized=True)
        ax_scatter.plot(xaxis, y2, 'dodgerblue')
        pnscatter = ax_scatter.scatter(Pn_x, Pn_y, c='darkred')
        p2scatter = ax_scatter.scatter(P2_x, P2_y, c='darkred')


        fig4.savefig(os.path.join(savefolder, "Phasor_2species_Confocal_{}.pdf".format(os.path.basename(mixedimage).split(".msr")[0])), transparent='True',
                     bbox_inches="tight",dpi=900)


        pnscatter.remove()
        p2scatter.remove()
        ax_scatter.lines[-1].remove()
        mixphasor.remove()


        lifetime=numpy.dstack(( fraction2, fraction1, fraction1*0))


        red=lifetime[:,:,1]+lifetime[:,:,2]
        green=lifetime[:,:,0]+lifetime[:,:,2]
        blue=lifetime[:,:,0]+lifetime[:,:,1]

        lifetime_rgb =numpy.dstack((red,green,blue))
        #lifetime_rgb=numpy.where(lifetime_rgb>1, 1,lifetime_rgb)
        print(lifetime_rgb.shape,lifetime_rgb.dtype)
        print(numpy.min(lifetime_rgb[:,:,0]),numpy.max(lifetime_rgb[:,:,0]))  
        print(numpy.min(lifetime_rgb[:,:,1]),numpy.max(lifetime_rgb[:,:,1]))  
        print(numpy.min(lifetime_rgb[:,:,2]),numpy.max(lifetime_rgb[:,:,2]))  

        plt.figure()
        plt.imshow(lifetime_rgb)
        plt.axis('off')

        plt.savefig(os.path.join(savefolder, os.path.basename(mixedimage).split(".msr")[0] + "_STED2species_lifetimergb.pdf"), transparent='True', bbox_inches="tight")
        filenameout = os.path.join(savefolder,
                                   os.path.basename(mixedimage).split(".msr")[0] + "_Conf2species_LineControls_lifetimergb.tiff")

        tifffile.imwrite(filenameout, lifetime_rgb.astype(numpy.float32))




        fig_im3, ax_im3 = plt.subplots(ncols=4,nrows=1,figsize=(12,8))
        ax_im3[0].axis('off')
        ax_im3[1].axis('off')
        ax_im3[2].axis('off')
        ax_im3[3].axis('off')
        #ax_im4.axis('off')
        #imsum_flat1 =ax_im.imshow(imsum_flat, cmap='turbo', vmin=numpy.min(t2), vmax=1)
        imsum_flat3 =ax_im3[1].imshow(imsum_flat_lin1, cmap='hot')
        cbar3 =fig_im3.colorbar(imsum_flat3, ax=ax_im3[1],fraction=0.05, pad=0.01)
        imsum_flat5 =ax_im3[2].imshow(imsum_flat_lin2, cmap='hot')
        cbar2 =fig_im3.colorbar(imsum_flat5, ax=ax_im3[2],fraction=0.05, pad=0.01)
        imsum_flat6 =ax_im3[0].imshow(imsum, cmap='hot')
        cbar1 =fig_im3.colorbar(imsum_flat6, ax=ax_im3[0],fraction=0.05, pad=0.01)
        overlayer = LifetimeOverlayer(fraction1, imsum/imsum.max(), cname='coolspring')
        lifetime_rgb, cmap = overlayer.get_overlay(
            lifetime_minmax=(0., 1),
            intensity_minmax=(0, 0.5) # inTensity saturated to get more bright regions
                    )
        imsum_flat5 =ax_im3[3].imshow(lifetime_rgb)
        cbar =fig_im3.colorbar(cmap, ax=ax_im3[3],fraction=0.05, pad=0.01)


        imagecomp=numpy.dstack((imsum_flat_lin2,imsum_flat_lin1))
        imagecomp=numpy.moveaxis(imagecomp,2,0)
        filenameout =  os.path.join(savefolder,os.path.basename(mixedimage).split(".msr")[0] + "_ConfCentroidsCircle_UnmixedComposite.tiff")
        imsave(file=filenameout, data=imagecomp.astype(numpy.uint16), composite=True, luts=("Cyan Hot","Magenta Hot"), pixelsize=(20E-3,20E-3))

        filenameout = os.path.join(savefolder,os.path.basename(mixedimage).split(".msr")[0] + "_ConfCentroidsCircle_MixedIntensity.tiff")
        print(filenameout)
        imsave(file=filenameout, data=imsum.astype(numpy.uint16), luts="Red Hot", pixelsize=(20E-3,20E-3))


        filenameout =  os.path.join(savefolder,os.path.basename(mixedimage).split(".msr")[0] + "_ConfCentroidsCircle_f1f2.tiff")
        imagecomp=numpy.dstack((fraction2,fraction1))
        imagecomp=numpy.moveaxis(imagecomp,2,0)
        tifffile.imwrite(filenameout, imagecomp)

        filenameout = os.path.join(savefolder,
                                   os.path.basename(mixedimage).split(".msr")[0] + "_Conf2species_LineControls_F1Overlay.tiff")

        tifffile.imwrite(filenameout, lifetime_rgb.astype(numpy.float32))




        fig_im3.savefig(os.path.join(savefolder,'Images__SeparateConf_CentroidsCircle'+os.path.basename(mixedimage).split(".msr")[0] +'.pdf'),transparent='True', bbox_inches="tight")

#plt.show()

