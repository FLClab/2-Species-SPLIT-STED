# -*- coding: utf-8 -*-
import os.path
from sys import path as path1;
dossier = os.path.expanduser("~/Documents/Github/2-Species-SPLIT-STED/Functions")
path1.append(dossier)
import math
import matplotlib.pyplot as plt

import numpy
import glob
import itertools
import seaborn
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import matplotlib.patches as mpatches
import pandas as pd
from Main_functions import (load_image,select_channel,line_equation, to_polar_coord, polar_to_cart, get_foreground)
from Phasor_functions import Median_Phasor
from Mono_fit import ExpFunMono_MLE
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression
import skimage
from lifetime import LifetimeOverlayer
from skimage import filters
import scipy
import decorr
from objectives import (Squirrel, Bleach)
import tifffile
from scipy.optimize import minimize
from tqdm import tqdm
matplotlib.rcParams['axes.linewidth'] = 0.8
# ------------------ Default Input variables ----------------
params_dict = {
    # Parameter in option in the matlab code
    #    "Tg" : 6, #% 'First frame to sum:'
    "Nb_to_sum": 100,  # The Tg infered from this variable override Tg
    "smooth_factor":0.2,  # % 'Smoothing factor:'
    "im_smooth_cycles": 0,  # % 'Smoothing cycles image:'
    "phasor_smooth_cycles": 1,  # % 'Smoothing cycles phasor:'
    "foreground_threshold": 10,
    "tau_exc": numpy.inf,  # % 'Tau_exc'
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
def Simulate2SpeciesSTED(STEDPOWER):

    f1 = os.path.join('T:', os.sep,'adeschenes','SimulationDataset_STEDFLIM','Cy3','Bassoon_CF594',"MediumAcq","*_{}percentSTED.msr".format(STEDPOWER))
    #f1 = os.path.join('U:', os.sep, 'adeschenes', "2023-12-21_FLIM_MediumAcq_Spectrin_Actin_Bassoon", "Bassoon_CF594_STEDPowerBleach_MediumAcq_1","*_{}percentSTED.msr".format(STEDPOWER))
    f2 = os.path.join('T:', os.sep,'adeschenes','SimulationDataset_STEDFLIM','Cy3','Homer_STORANGE',"MediumAcq","*_{}percentSTED.msr".format(STEDPOWER))
    f2 = os.path.join('T:', os.sep, 'adeschenes', 'SimulationDataset_STEDFLIM', 'Cy3', "PSD95_STORANGE","Mini20",
                      "*_{}percentSTED.msr".format(STEDPOWER))
    f1 = os.path.join('T:', os.sep, 'adeschenes', 'SimulationDataset_STEDFLIM', 'Cy3', "rabBassoon_CF594","Mini20",
                      "*_{}percentSTED.msr".format(STEDPOWER))

    filenames = [f1,f2]
    filenamescontrol = [f1,f2]
    labels=["Bassoon_CF594",'PSD95_STORANGE',"Mixture"]
    #keys=['CONF640 {10}','CONF640 {0}']
    keyscontrols = ['STED 561 {11}','STED 561 {11}']
    keysmixed = ['STED 561 {11}','STED 561 {11}']
    #keyscontrols=[1,1]
    #keysmixed=[1,1]
    #keyscontrols = [ 'Conf_635P {2}','Conf_635P {2}']
    #keysmixed = [ 'STED_635P {2}', 'STED_635P {2}']
    
    msrfiles = []
    #plt.style.use('dark_background')
    
    colors=['cyan','magenta','springgreen']
    #colors=['magenta']
    
    
    
    
    
    #savefolder=str(input("Name of Output folder: "))
    savefolder = "Simulation_Cy3_{}Percent_MLEinput".format(STEDPOWER)
    savefolder = os.path.join(os.path.expanduser("~/Desktop"), savefolder)
    os.makedirs(savefolder,exist_ok=True)

 
    msrfiles=[]
    images=[glob.glob(filename)for filename in filenames]
    number = [len(glob.glob(filename)) for filename in filenames]
    print('There are ',number, 'Images in these folders')
    pairs = list(itertools.product(images[0], images[1]))
    print(len(pairs))
    Overall_data = pd.DataFrame(
        columns=["Power",'image1', 'image2', 'resolution1', 'resolution2'])
    
    fig4,ax_scatter = plt.subplots(figsize=(2,2))

    edge = numpy.linspace(start=0, stop=15, num=200)
    theta = numpy.linspace(0, numpy.pi, 100)
    r = 0.5
    x1 = r * numpy.cos(theta) + 0.5
    x2 = r * numpy.sin(theta)

    ax_scatter.plot(x1, x2, color="k", ls="--",linewidth=0.8)
    ax_scatter.set_xlim(0, 1)
    ax_scatter.set_ylim(0, 1)
    ax_scatter.set_xlabel('g')
    ax_scatter.set_ylabel('s')

    for Pair_id, (a, b) in enumerate(pairs):
        ov_data = [STEDPOWER,a, b]
        msrfiles = [a, b]
        print("***********************************************************")
        print("Working on pair number ",Pair_id," out of ",len(pairs))
        print("***********************************************************")
        print(msrfiles)
        croplist=[]
        Imagelist = []
        Combolist = []
        ComboMasklist = []
        Masklist = []
        ControlImagesList = []
        seuils = []
        for i, msr in enumerate(msrfiles):
            print("i",i)
            imagemsr=load_image(msr)
            #print(imagemsr.keys())
            image1 = select_channel(imagemsr, keysmixed[i])
            imagec1 = select_channel(imagemsr, keyscontrols[i])
            #image1 = imagemsr[keysmixed[i]]
            #imagec1 = imagemsr[keyscontrols[i]]
            #res_control = decorr.calculate(numpy.sum(imagec1[:,:,10:],axis=2))
            res_mix = decorr.calculate(numpy.sum(image1[:, :, 10:111],axis=2))
            if math.isinf(res_mix):
                res_mix=10
            #print("res_control",res_control*20,"res_mix ",res_mix*20 )
            ov_data.append(res_mix * 20)
            #ControlImagesList.append(imagec1)
            Imagelist.append(image1)
            print(image1.shape)
            imsum=numpy.sum(image1[:,:,10:111],axis=2)
            #seuil = get_foreground(image1)
            seuil=10
            seuils.append(seuil)
            mask=imsum>seuil
            mask=scipy.ndimage.binary_fill_holes(mask)
            Masklist.append(mask)
    
        print("masklist",len(Masklist),Masklist[0].shape)
        Combo=numpy.zeros(numpy.max([Imagelist[0].shape,Imagelist[1].shape],axis=0))
        Comboo = numpy.zeros(numpy.max([Imagelist[0].shape, Imagelist[1].shape], axis=0))
        #Combo=numpy.zeros((300,300,250))
        Combomask1 = numpy.zeros(Combo.shape[0:2])
        Combomask2 = numpy.zeros(Combo.shape[0:2])
        Combosingle1 = numpy.zeros(Combo.shape[0:2])
        Combosingle2= numpy.zeros(Combo.shape[0:2])
        print('Combo',Combo.shape)
        print('Combomask', Combomask1.shape)
    
        minx=0
        maxx=Imagelist[0].shape[0]
        miny=0
        maxy=Imagelist[0].shape[1]

        Combo[minx:maxx,miny:maxy ,:]+=Imagelist[0][:,:, :]
        Combomask1[minx:maxx,miny:maxy]+=Masklist[0]
        Combolist.append(Combo.copy())
        Combosingle1[minx:maxx, miny:maxy] += numpy.sum(Imagelist[0][:,:,10:111],axis=2)
        minx = 0
        maxx = Imagelist[1].shape[0]
        miny = 0
        maxy = Imagelist[1].shape[1]
        print(minx, miny, maxx, maxy)
    
        Combo[minx:maxx, miny:maxy, :] += Imagelist[1][:, :, :]
        Comboo[minx:maxx, miny:maxy, :] += Imagelist[1][:, :, :]
        Combomask2[minx:maxx, miny:maxy] += Masklist[1]
    
        Combosingle2[minx:maxx, miny:maxy] += numpy.sum(Imagelist[1][:,:,10:111], axis=2)
        ComboMasklist = [Combomask1, Combomask2]
        Combolist.append(Comboo)
        Combolist.append(Combo)

    
        Imagelist.append(Combo)
        ControlImagesList.append(Combo)
    
        CoM_x, CoM_y = [], []
        d_melange = pd.DataFrame(columns=['g', 's'])
        for i, image1 in enumerate(Combolist):
            df = pd.DataFrame(columns=['x', 'y'])
            dg = pd.DataFrame(columns=['g', 's'])
    
            print(image1.shape)
            imsum = image1[:,:,10:111].sum(axis=2)
            imsum = imsum.astype('int16')
    
            #seuil = get_foreground(imsum)
            #seuil=min(seuils)
            seuil=10
            print("Caclulation for an image of shape", image1.shape, "...")
    
            params_dict["foreground_threshold"] = seuil
            params_dict["Nb_to_sum"] = image1.shape[2]
            print("foreground_threshold=", params_dict["foreground_threshold"])
            #x,y, original_idxes,Images,Images_Filtered=DTCWT_Phasor(image1, 0, nlevels=10, neighborhood=50)
            x, y, g_smoothed, s_smoothed, original_idxes = Median_Phasor(image1, params_dict, **params_dict, show_plots=False)
            #x = x[imsum > params_dict["foreground_threshold"]]
            #y = y[imsum > params_dict["foreground_threshold"]]
            df['x'] = x.flatten()
            df['y'] = y.flatten()
            m, phi = to_polar_coord(df['x'], df['y'])
            g, s = polar_to_cart(m, phi)
            dg['g'], dg['s'] = g, s
    
            kmeans = KMeans(n_clusters=1, init='k-means++', random_state=42)
            y_kmeans = kmeans.fit_predict(dg)
            CoM_x.extend(kmeans.cluster_centers_[:, 0][:].tolist())
            CoM_y.extend(kmeans.cluster_centers_[:, 1][:].tolist())
            mixphasor = ax_scatter.scatter(g, s, c=colors[i],s=0.5,alpha=0.1,rasterized=True)

            fig4.savefig(os.path.join(savefolder,"Phasor_{}_{}.pdf".format(Pair_id,labels[i])),transparent='True', bbox_inches="tight",dpi=900)

            mixphasor.remove()


            seuil= 10
            mlifetime = numpy.empty((imsum.shape))

            for iy, ix in tqdm(numpy.ndindex(imsum.shape)):
                y = image1[iy, ix]
                if y.sum() < seuil :
                    mlifetime[iy, ix] = 0
                else :
                    maxy = numpy.max(y) 
                    indice = numpy.argmax(y)
                    y = y[indice:]
                    y= y / y.sum()
                    absci = numpy.linspace(0,y.shape[0]-1, num =y.shape[0])*0.08
        #   Computes lifetime for each pixel

                    tau = 2.5
                    nb_photons = 100
                    bounds = [(0,0.0001),(0.1,5),(0,900000)] # (min, max): Background, tau, amp
                    w = minimize(ExpFunMono_MLE, 
                                x0 = [0,tau,nb_photons], 
                                args = [absci, y], 
                                bounds = bounds)
                    bg, tau, amp =  w["x"]
                    mlifetime[iy, ix] = w['x'][1]

            fig1, ax1 = plt.subplots()
            imgplot1 = ax1.imshow(mlifetime, cmap='jet')
            ax1.axis('off')
            cbar =fig1.colorbar(imgplot1)
            cbar.set_label("Lifetime [ns]")
            #plt.show()

        # -----------------------------------------------------------
        #    Lifetime modulated in intensity


            overlayer = LifetimeOverlayer(mlifetime, imsum/imsum.max(), cname='jet')
            lifetime_rgb, cmap = overlayer.get_overlay(
                lifetime_minmax=(1,4),
                intensity_minmax=(0, 0.3)
                        )
                    
            """ If error, watch out conversion in lifetime.py """

            fig2, ax2 = plt.subplots()
            ax2.axis('off')
            img = ax2.imshow(lifetime_rgb)
            cbar = fig2.colorbar(cmap, ax=ax2)
            cbar.set_label("temps de vie [ns]")
            filenameout = os.path.join(savefolder,"MLE_Lifetime_{}_{}.tiff".format(Pair_id,labels[i]))
            tifffile.imwrite(filenameout, mlifetime)
            filenameout =os.path.join(savefolder,"Intensity_{}_{}.tiff".format(Pair_id,labels[i]))
            tifffile.imwrite(filenameout, imsum.astype(numpy.uint16))
            fig1.savefig(os.path.join(savefolder,"MLELifetime_{}_{}.pdf".format(Pair_id,labels[i])), transparent='True', bbox_inches="tight")
            fig2.savefig(os.path.join(savefolder,"MLELifetime_IntensityComposite_{}_{}.pdf".format(Pair_id,labels[i])), transparent='True', bbox_inches="tight")

            plt.close(fig1)
            plt.close(fig2)
    
        GroundTruth_Fraction=[(Combosingle1/(Combosingle2+Combosingle1)),(Combosingle2 / (Combosingle2 + Combosingle1))]
        imagecomp=numpy.dstack((GroundTruth_Fraction[0],GroundTruth_Fraction[1],Combosingle1,Combosingle2,imsum))
        imagecomp=numpy.moveaxis(imagecomp,2,0)
        filenameout =os.path.join(savefolder, "{}_GroundTruth.tiff".format(Pair_id))
        print(filenameout)
        tifffile.imwrite(filenameout, imagecomp)
        with open(os.path.join(savefolder,'legend.txt'),'a') as data:
            data.write("{}\t{}\t{}\n".format(Pair_id,a,b))




    return numpy.array([Overall_data["resolution1"],Overall_data["resolution2"]])
#["_30",[0,4]],["30",[0,0]],["40",[0,0]]

#[[5,[0,0]],[10,[0,0]],[15,[0,0]],[20,[0,0]],[30,[0,0]],[40,[0,0]]]
Powerslist=[10,20,30,40]
Powerslist=[20]
#Powerslist=[[30,[0,0]],[40,[0,0]]]
globalcumstats=[]
globalcumstatsmean=[]
globalcumstatsstd=[]
row=0

for power in Powerslist:
    stats=Simulate2SpeciesSTED(power)
    print(stats.shape)

    plt.show()