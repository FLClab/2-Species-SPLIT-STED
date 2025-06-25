# -*- coding: utf-8 -*-
import os.path
from sys import path as path1;

Functionspath=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Functions")
path1.append(Functionspath)
import easygui
from Main_functions import (line_equation, to_polar_coord, polar_to_cart, load_image,select_channel, get_foreground)
from Phasor_functions import Median_Phasor,DTCWT_Phasor,unmix2species
from tiffwrapper import imsave,LifetimeOverlayer
from objectives import (Squirrel, Bleach)
import math
import matplotlib.pyplot as plt
import numpy
import glob
import itertools
import seaborn
import matplotlib
import matplotlib.patches as mpatches
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression
import skimage
from skimage import filters
import scipy
import decorr
import tifffile
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
def Simulate2SpeciesSTED(STEDPOWER,NUMIM):

    f1=easygui.diropenbox(default=os.path.expanduser("~Desktop"),title="Select folder containing control images for First fluorophore")
    f2=easygui.diropenbox(default=os.path.expanduser("~Desktop"),title="Select folder containing control images for Second fluorophore")

    f2= os.path.join(f2,"*_{}percentSTED.msr".format(STEDPOWER))
    f1= os.path.join(f1,"*_{}percentSTED.msr".format(STEDPOWER))
 
    filenames = [f1,f2]
    filenamescontrol = [f1,f2]
    labels=["Bassoon_CF594",'PSD95_STORANGE',"Mixture"]
  
    keyscontrols = ['STED 561 {11}','STED 561 {11}']
    keysmixed = ['STED 561 {11}','STED 561 {11}']

    #keyscontrols = ['STED_635P {2}', 'STED_635P {2}']
    #keysmixed = [ 'STED_635P {2}', 'STED_635P {2}']
    
    msrfiles = []
    #plt.style.use('dark_background')
    
    colors=['royalblue','orangered','springgreen']
    #colors=['magenta']
    
    
    
    
    
    #savefolder=str(input("Name of Output folder: "))
    savefolder = "Simulation_Cy3_{}Percent_2Species_PSD95Bassoon".format(STEDPOWER)
    savefolder = os.path.join(os.path.expanduser("~/Desktop"), savefolder)
    os.makedirs(savefolder,exist_ok=True)

    for k, filename in enumerate(filenamescontrol):
        print(labels[k])
        #path = os.path.join(filename, '*.msr')
        images = glob.glob(filename)
        print('There are ', len(images), 'Images in this folder')
        for imagei in images:
            print(os.path.basename(imagei))
        #numim = int(input('Fichier msr a extraire (1er=0): '))
        numim = NUMIM[k]
        image = images[numim]
        msrfiles.append(image)
    print(msrfiles)
    
    CoM_x, CoM_y = [], []
    
    fig4,ax_scatter = plt.subplots(figsize=(3,3))

    edge = numpy.linspace(start=0, stop=15, num=200)
    theta = numpy.linspace(0, numpy.pi, 100)
    r = 0.5
    x1 = r * numpy.cos(theta) + 0.5
    x2 = r * numpy.sin(theta)
    ax_scatter.plot(x1, x2, color="k", ls="--",linewidth=0.8)
    
    with open(os.path.join(savefolder,'legend.txt'),'w') as data: 
        data.write("Controls\n")
    scatterlist=[]

    for i, msr in enumerate(msrfiles) : 
        df = pd.DataFrame(columns=['x','y'])
        dg = pd.DataFrame(columns=['g', 's'])
        
        with open(os.path.join(savefolder,'legend.txt'),'a') as data: 
            data.write("{}\t{}\t{}\n".format(labels[i],keyscontrols[i],msr))
        imagemsr = load_image(msr)
        image1= select_channel(imagemsr, keyscontrols[i])
        #image1 = imagemsr[keyscontrols[i]]
        imsum = image1[:,:,10:111].sum(axis=2)
        imsum = imsum.astype('int16')
    
        seuil = get_foreground(imsum)
        print("Caclulation for an image of shape", image1.shape, "...")
        params_dict["foreground_threshold"] = seuil
        params_dict["Nb_to_sum"] = image1.shape[2]
        print("foreground_threshold=", params_dict["foreground_threshold"])
    
        x, y, g_smoothed, s_smoothed, original_idxes = Median_Phasor(image1, params_dict, **params_dict,
                                                                     show_plots=False)
        df['x'] = x.flatten()
        df['y'] = y.flatten()
        m, phi = to_polar_coord(df['x'], df['y'])
        g, s = polar_to_cart(m, phi)
        dg['g'], dg['s'] = g, s
        kmeans = KMeans(n_clusters=1, init='k-means++', random_state=42)
        y_kmeans = kmeans.fit_predict(dg)
        CoM_x.extend(kmeans.cluster_centers_[:, 0][:].tolist())
        CoM_y.extend(kmeans.cluster_centers_[:, 1][:].tolist())
        a=ax_scatter.scatter(g, s, s=1, c=colors[i], alpha=0.10)
        scatterlist.append(a)
        #print('DENSITY ESTIMATION IS HARD WORK, BE PATIENT PLEASE')
        #seaborn.kdeplot(x=g,y=s,ax=ax_scatter,color=colorsc[i],levels=[0.2,0.4,0.6,0.8,1.0],linewidths= 1.5,label=labels[i])


    ax_scatter.set_xlim(0, 1)
    ax_scatter.set_ylim(0, 1)

    ax_scatter.set_xlabel('g')
    ax_scatter.set_ylabel('s')
    #fig4.savefig(os.path.join(savefolder, "Phasor_2species_ControlsOnly.png"), transparent='True', bbox_inches="tight")
    xaxis = numpy.linspace(0, 1, 100)

    
    ##Calcul de Pn
    Pn_x, Pn_y = CoM_x[0], CoM_y[0]
    P_n = numpy.array([Pn_x, Pn_y])
    P2_x, P2_y = CoM_x[1], CoM_y[1]
    p2 = numpy.array([P2_x, P2_y])
    ## Droite entre Pn - p2
    x2, y2 = [Pn_x, P2_x], [Pn_y, P2_y]
    m2, c2 = line_equation(x2, y2)
    y2 = m2 * xaxis + c2
    

    #lines = [mpatches.Patch(color=colors[j], label=labels[j]) for j in range(len(labels))]
    #ax_scatter.legend(handles=lines, prop={'size': 20})
    t = [scatter.remove() for scatter in scatterlist]
    pnscatter = ax_scatter.scatter(Pn_x, Pn_y, s=50, c='gold')
    p2scatter = ax_scatter.scatter(P2_x, P2_y, s=50, c='gold')
    p2pnline = ax_scatter.plot([Pn_x, P2_x], [Pn_y, P2_y], c='dodgerblue')

    fig4.savefig(os.path.join(savefolder, "Phasor_2species_ControlsOnly.pdf"), transparent='True', bbox_inches="tight")
    #ax_scatter.get_legend().remove()
    pnscatter.remove()
    p2scatter.remove()
    ax_scatter.lines[-1].remove()


    #fig4.savefig(os.path.join(savefolder, "Phasor_2species_ControlsOnly.png"), transparent='True', bbox_inches="tight")
    msrfiles=[]
    images=[glob.glob(filename)for filename in filenames]
    number = [len(glob.glob(filename)) for filename in filenames]
    print('There are ',number, 'Images in these folders')
    pairs = list(itertools.product(images[0], images[1]))
    print(len(pairs))
    Overall_data = pd.DataFrame(
        columns=["Power",'image1', 'image2', 'resolution1', 'resolution2', 'TrueFraction1', 'PredictedFraction1', 'PixelError1',
                 'PixelIntensity1', 'fit_intercept1', 'fit_slope1', 'fit_correlation1',  'TrueFraction2', 'PredictedFraction2', 'PixelError2',
                 'PixelIntensity2', 'fit_intercept2', 'fit_slope2', 'fit_correlation2',"res_fraction1", "res_fraction2","squirrel_f1","squirrel_f2",
                 "squirrelsmooth_f1","squirrelsmooth_f2"])
    for Pair_id, (a, b) in enumerate(pairs):
        ov_data = [STEDPOWER,a, b]
        msrfiles = [a, b]
        print("***********************************************************")
        print("Working on pair number ",Pair_id," out of ",len(pairs))
        print("***********************************************************")
        print(msrfiles)
        croplist=[]
        Imagelist = []
        data = pd.DataFrame(columns=['area', 'bbox', 'mean_intensity', 'coords', 'centroid', 'channel'])
        ComboMasklist = []
        Masklist = []
        Propslist = []
        Cropslist = []
        CropImageList = []
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
            seuil=3
            seuils.append(seuil)
            mask=imsum>seuil
            mask=scipy.ndimage.binary_fill_holes(mask)
            Masklist.append(mask)
    
        print("masklist",len(Masklist),Masklist[0].shape)
        Combo=numpy.zeros(numpy.max([Imagelist[0].shape,Imagelist[1].shape],axis=0))
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
        print(minx, miny, maxx, maxy)
    
    
        #CropImageList.append(numpy.sum(Imagelist[0][bbox[0]:bbox[2],bbox[1]:bbox[3], :]*slice3d,axis=2))
        Combo[minx:maxx,miny:maxy ,:]+=Imagelist[0][:,:, :]
        Combomask1[minx:maxx,miny:maxy]+=Masklist[0]
    
        Combosingle1[minx:maxx, miny:maxy] += numpy.sum(Imagelist[0][:,:,10:111],axis=2)
        minx = 0
        maxx = Imagelist[1].shape[0]
        miny = 0
        maxy = Imagelist[1].shape[1]
        print(minx, miny, maxx, maxy)
    
        Combo[minx:maxx, miny:maxy, :] += Imagelist[1][:, :, :]
        Combomask2[minx:maxx, miny:maxy] += Masklist[1]
    
        Combosingle2[minx:maxx, miny:maxy] += numpy.sum(Imagelist[1][:,:,10:111], axis=2)
        ComboMasklist = [Combomask1, Combomask2]
        ComboSinglelist=[Combosingle1,Combosingle2]
    
        #Combo = filters.gaussian(Combo, sigma)
        Imagelist.append(Combo)
        ControlImagesList.append(Combo)
    
        CoM_x, CoM_y = [], []
        d_melange = pd.DataFrame(columns=['g', 's'])
        for i, image1 in enumerate(ControlImagesList):
            df = pd.DataFrame(columns=['x', 'y'])
            dg = pd.DataFrame(columns=['g', 's'])
    
            print(image1.shape)
            imsum = image1[:,:,10:111].sum(axis=2)
            imsum = imsum.astype('int16')
    
            #seuil = get_foreground(imsum)
            seuil=min(seuils)
            #seuil=10
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

    
        p3 = dg[['g', 's']].to_numpy() #phaseur qui sera projeté
        imsum_flat_lin1, imsum_flat_lin2, Solve = unmix2species(p3, original_idxes, Imagelist[2], P_n, p2)
    
        fraction2 = imsum_flat_lin1.copy()
        fraction1 = imsum_flat_lin2.copy()
    
    
        GroundTruth_Fraction=[(Combosingle1/(Combosingle2+Combosingle1)),(Combosingle2 / (Combosingle2 + Combosingle1))]
        Predicted_Fraction=[fraction2.copy(),fraction1.copy()]

        #lines = [mpatches.Patch(color=colors[j], label=labels[j]) for j in range(len(labels))]
        #ax_scatter.legend(handles=lines, prop={'size': 20})

        mixphasor = ax_scatter.scatter(g, s, s=2,c=Solve[1,:],cmap="cool",rasterized=True,label="Mixture")
        lineplot=ax_scatter.plot(xaxis, y2, 'dodgerblue')
        pnplot=ax_scatter.scatter(Pn_x, Pn_y, s=50,c='gold')
        p2plot=ax_scatter.scatter(P2_x, P2_y, s=50,c='gold')
        fig4.savefig(os.path.join(savefolder,"Phasor_2species_{}.pdf".format(Pair_id)),transparent='True', bbox_inches="tight",dpi=900)
        #ax_scatter.get_legend().remove()
        pnplot.remove()
        p2plot.remove()
        ax_scatter.lines[-1].remove()
        #fig4.savefig(os.path.join(savefolder, "Phasor_2species_{}.png".format(Pair_id)), transparent='True',
        #                bbox_inches="tight",dpi=900)
        mixphasor.remove()
    
        FilteredMaskList=[]
        for p,maski in enumerate(ComboMasklist):
            TruefractionPixels=[]
            PredfractionPixels=[]
            IntensityPixels=[]
            Truefraction=[]
            Predfraction=[]
            labelsim = skimage.measure.label(numpy.logical_or(ComboMasklist[0],ComboMasklist[1]))
            props = skimage.measure.regionprops(labelsim, intensity_image=Combo[:, :, 10:111].sum(axis=2))
            mask = numpy.zeros(Combo.shape[0:2], dtype=bool)
            print(ComboSinglelist[p].shape,numpy.min(maski*ComboSinglelist[p]),numpy.max(maski*ComboSinglelist[p]))
            for j, region in enumerate(props):
                #if (region.area > 10) and (0 not in region.bbox) and (Imagelist[p].shape[0] not in region.bbox) and (Imagelist[p].shape[1] not in region.bbox):
                if (region.area > 10) :
                    data.loc[data.shape[0] + 1] = {'area': region.area, 'bbox': region.bbox,
                                               'mean_intensity': region.mean_intensity, 'coords': region.image,
                                               'centroid': region.centroid, 'channel': p}
                    # print(region.image)
                    croplist.append(region.intensity_image)
                    mask[region.bbox[0]:region.bbox[2], region.bbox[1]:region.bbox[3]] += region.image
                    truefraction=GroundTruth_Fraction[p][region.bbox[0]:region.bbox[2], region.bbox[1]:region.bbox[3]]*region.image
                    predfraction = Predicted_Fraction[p][region.bbox[0]:region.bbox[2],
                                   region.bbox[1]:region.bbox[3]] * region.image
                    #print(numpy.count_nonzero(numpy.isnan(truefraction)))
                    TruefractionPixels.extend( truefraction[numpy.isfinite(truefraction) ].flatten())
                    #print(numpy.nanmean(truefraction))
                    Truefraction.append(numpy.nanmean(truefraction))
    
                    Predfraction.append(numpy.mean(predfraction))
                    PredfractionPixels.extend(predfraction[numpy.isfinite(truefraction)].flatten())
                    IntensityPixels.extend((numpy.sum(Combo[:,:,10:111],axis=2)[region.bbox[0]:region.bbox[2], region.bbox[1]:region.bbox[3]]*region.image)[numpy.isfinite(truefraction) ].flatten())
            FilteredMaskList.append(mask)
        
            TruefractionPixels = numpy.array(TruefractionPixels)
            PredfractionPixels=numpy.array(PredfractionPixels)
            ErrorPixels=numpy.abs(PredfractionPixels-TruefractionPixels)
            IntensityPixels = numpy.array(IntensityPixels)
    
            print(data['channel'].unique())
            model = LinearRegression().fit(TruefractionPixels.reshape((-1, 1)),PredfractionPixels)
            r_sq = model.score(TruefractionPixels.reshape((-1, 1)),PredfractionPixels)
            x_pred=numpy.linspace(0, 1, 10)
            y_pred = model.intercept_ + model.coef_ * x_pred
    
            ov_data.extend(
                [TruefractionPixels[numpy.nonzero(IntensityPixels)], PredfractionPixels[numpy.nonzero(IntensityPixels)],
                ErrorPixels[numpy.nonzero(IntensityPixels)], IntensityPixels[numpy.nonzero(IntensityPixels)], model.intercept_,
                model.coef_[0], r_sq])
        
    
        Combomask_nofilter=numpy.logical_or(ComboMasklist[0],ComboMasklist[1])
        Combomask = numpy.logical_or(FilteredMaskList[0],FilteredMaskList[1])
        imsum = Combo[:, :, 10:111].sum(axis=2)
        fraction1 *= imsum
        fraction2 *= imsum
        res_fraction1 = decorr.calculate(fraction1.astype(numpy.uint16) )
        res_fraction2 = decorr.calculate(fraction2.astype(numpy.uint16) )
        if math.isinf(res_fraction2):
            res_fraction2=10
        if math.isinf(res_fraction1):
            res_fraction1=10
        print("res_fraction1", res_fraction1 * 20, "res_fraction2", res_fraction2 * 20)

        squirrel_f1= Squirrel(method="L-BFGS-B", normalize=True).evaluate([fraction2], Combosingle1*Combomask1,
                                                                                      Combosingle1*Combomask1,
                                                                                      Combomask1,Combomask1)
        squirrel_f2= Squirrel(method="L-BFGS-B", normalize=True).evaluate([fraction1], Combosingle2*Combomask2,
                                                                                      Combosingle2*Combomask2,
                                                                                      Combomask2,Combomask2)
        squirrelmap1,squirrelsmoothf1 = Squirrel(method="L-BFGS-B", normalize=True).return_map([fraction2], Combosingle1 * Combomask1,
                                                                           Combosingle1 * Combomask1,
                                                                           Combomask1, Combomask1)
        squirrelmap2 ,squirrelsmoothf2= Squirrel(method="L-BFGS-B", normalize=True).return_map([fraction1], Combosingle2 * Combomask2,
                                                                           Combosingle2 * Combomask2,
                                                                           Combomask2, Combomask2)


        ov_data.extend([res_fraction1* 20, res_fraction2* 20,squirrel_f1,squirrel_f2,squirrelsmoothf1[2],squirrelsmoothf2[2]])
        Overall_data.loc[Pair_id] = ov_data

        imagecomp=numpy.dstack((squirrelmap1,squirrelmap2))
        imagecomp=numpy.moveaxis(imagecomp,2,0)
        filenameout = os.path.join(savefolder,"{}_SquirrelMaps.tiff".format(Pair_id))
        tifffile.imwrite(filenameout, imagecomp)
        
        imagecomp=numpy.dstack((GroundTruth_Fraction[0],GroundTruth_Fraction[1],Combosingle1,Combosingle2,imsum))
        imagecomp=numpy.moveaxis(imagecomp,2,0)
        filenameout =os.path.join(savefolder, "{}_GroundTruth.tiff".format(Pair_id))
        print(filenameout)
        tifffile.imwrite(filenameout, imagecomp)

        imagecomp=numpy.dstack((ComboMasklist[0],ComboMasklist[1],Combomask_nofilter,FilteredMaskList[0],FilteredMaskList[1],Combomask))
        imagecomp=numpy.moveaxis(imagecomp,2,0)
        filenameout = os.path.join(savefolder,"{}_Masks.tiff".format(Pair_id))
        tifffile.imwrite(filenameout, imagecomp)
    
        imagecomp=numpy.dstack((Predicted_Fraction[0],Predicted_Fraction[1],fraction1,fraction2))
        imagecomp=numpy.moveaxis(imagecomp,2,0)
        filenameout = os.path.join(savefolder,"{}_Predictions.tiff".format(Pair_id))
        tifffile.imwrite(filenameout, imagecomp)
    
        imagecomp=numpy.dstack((Combosingle1,Combosingle2))
        imagecomp=numpy.moveaxis(imagecomp,2,0)
        filenameout =os.path.join(savefolder, "{}_GroundTruth_Composite.tiff".format(Pair_id))
        imsave(file=filenameout, data=imagecomp.astype(numpy.uint16), composite=True, luts=("Cyan Hot","Magenta Hot"), pixelsize=(20E-3,20E-3))
    
        imagecomp=numpy.dstack((fraction2,fraction1))
        imagecomp=numpy.moveaxis(imagecomp,2,0)
        filenameout =os.path.join(savefolder, "{}_Predicted_Composite.tiff".format(Pair_id))
        imsave(file=filenameout, data=imagecomp.astype(numpy.uint16), composite=True, luts=("Cyan Hot","Magenta Hot"), pixelsize=(20E-3,20E-3))
       
        with open(os.path.join(savefolder,'legend.txt'),'a') as data: 
            data.write("{}\t{}\t{}\n".format(Pair_id,a,b))
    

            
        Predicted_Fraction[1][numpy.isnan(GroundTruth_Fraction[1])]=numpy.nan
        Predicted_Fraction[0][numpy.isnan(GroundTruth_Fraction[0])]=numpy.nan  
        fig_im, ax_im=plt.subplots(ncols=4, nrows=3, figsize=(12, 8), sharex=True, sharey=True)
        ax_im[0,0].axis('off')
        ax_im[0,1].axis('off')
        ax_im[0,2].axis('off')
        ax_im[0,3].axis('off')
        ax_im[1,0].axis('off')
        ax_im[1,1].axis('off')
        ax_im[1,2].axis('off')
        ax_im[1,3].axis('off')
        ax_im[2,0].axis("off")
        ax_im[2, 1].axis("off")
        ax_im[2,2].axis("off")
        ax_im[2, 3].axis("off")
    
        #
        ax_im[0,0].set_title('Initial image 1', fontsize=16)
        ax_im[0,1].set_title('Initial image 2', fontsize=16)
        ax_im[0,2].set_title('Ground Truth Mask 1', fontsize=16)
        ax_im[0,3].set_title('Ground Truth Mask 2', fontsize=16)
        ax_im[1,0].set_title('Fraction1', fontsize=16)
        ax_im[1,1].set_title('Fraction2', fontsize=16)
        ax_im[1,2].set_title('Predicted Fraction 1', fontsize=16)
        ax_im[1,3].set_title('Predicted Fraction 2', fontsize=16)
        ax_im[2,0].set_title('Error Fraction 1', fontsize=16)
        ax_im[2,1].set_title('Error Fraction 2', fontsize=16)
    
        imdisp0 = ax_im[0,0].imshow(Combosingle1, cmap='hot')
        imdisp1 = ax_im[0,1].imshow(Combosingle2, cmap='hot')
        cbar = fig_im.colorbar(imdisp0, ax=ax_im[0,0], fraction=0.05, pad=0.01)
        cbar = fig_im.colorbar(imdisp1, ax=ax_im[0,1], fraction=0.05, pad=0.01)
    
        imdisp0 = ax_im[1,0].imshow(fraction2, cmap='hot')
        imdisp1 = ax_im[1,1].imshow(fraction1, cmap='hot')
        cbar = fig_im.colorbar(imdisp0, ax=ax_im[1,0], fraction=0.05, pad=0.01)
        cbar = fig_im.colorbar(imdisp1, ax=ax_im[1,1], fraction=0.05, pad=0.01)
        imsum_flat3 = ax_im[0,3].imshow(GroundTruth_Fraction[1] * Combomask, cmap='rainbow')
        cbar3 = fig_im.colorbar(imsum_flat3, ax=ax_im[0,3], fraction=0.05, pad=0.01)
        imsum_flat3 = ax_im[1,2].imshow(Predicted_Fraction[0] * Combomask, cmap='rainbow')
        cbar4 = fig_im.colorbar(imsum_flat3, ax=ax_im[1,2], fraction=0.05, pad=0.01)
        imsum_flat3 = ax_im[0,2].imshow(GroundTruth_Fraction[0] * Combomask, cmap='rainbow')
        cbar1 = fig_im.colorbar(imsum_flat3, ax=ax_im[0,2], fraction=0.05, pad=0.01)
        imsum_flat5 = ax_im[1,3].imshow(Predicted_Fraction[1] *Combomask, cmap='rainbow')
        cbar2 = fig_im.colorbar(imsum_flat5, ax=ax_im[1,3], fraction=0.05, pad=0.01)
        imerr2 = ax_im[2,0].imshow(numpy.abs(Predicted_Fraction[0]-GroundTruth_Fraction[0])* Combomask ,cmap='Reds')
        cbarerr2 = fig_im.colorbar(imerr2 , ax=ax_im[2,0], fraction=0.05, pad=0.01)
        imerr1=ax_im[2,1].imshow(numpy.abs(Predicted_Fraction[1]-GroundTruth_Fraction[1]) *Combomask,cmap='Reds')
        cbarerr1 = fig_im.colorbar(imerr1, ax=ax_im[2,1], fraction=0.05, pad=0.01)
        imerr2 = ax_im[2, 2].imshow(GroundTruth_Fraction[0] ,cmap='Reds')
        cbarerr2 = fig_im.colorbar(imerr2 , ax=ax_im[2,2], fraction=0.05, pad=0.01)
        imerr1=ax_im[2,3].imshow(GroundTruth_Fraction[1],cmap='Reds')
        cbarerr1 = fig_im.colorbar(imerr1, ax=ax_im[2,3], fraction=0.05, pad=0.01)
        fig_im.savefig(os.path.join(savefolder,"Images_2species_{}.pdf".format(Pair_id)),transparent='True', bbox_inches="tight")
        plt.close(fig_im)
    


        overlayer = LifetimeOverlayer(fraction1, imsum/imsum.max(), cname='CET-I1')
        lifetime_rgb, cmap = overlayer.get_overlay(
             lifetime_minmax=(0., 1),
             intensity_minmax=(0, 0.5)  # inTensity saturated to get more bright regions
         )
        filenameout = os.path.join(savefolder,"{}_OverlayF1.tiff".format(Pair_id))
        tifffile.imwrite(filenameout, lifetime_rgb.astype(numpy.float32))    
        # overlayer = LifetimeOverlayer(fraction1, imsum/imsum.max(), cname='CET-D11')
        # lifetime_rgb, cmap = overlayer.get_overlay(
        #     lifetime_minmax=(0., 1),
        #     intensity_minmax=(0, 1) # inTensity saturated to get more bright regions
        #             )
        # imsum_flat5 =ax_im[1,4].imshow(lifetime_rgb)
        # cbar6 =fig_im.colorbar(cmap, ax=ax_im[1,4],fraction=0.05, pad=0.01)
    
    TrueFractionPixels=numpy.concatenate(Overall_data["TrueFraction1"].to_numpy())
    PredfractionPixels = numpy.concatenate(Overall_data["PredictedFraction1"].to_numpy())
    IntensityPixels= numpy.concatenate(Overall_data["PixelIntensity1"].to_numpy())
    ErrorPixels=numpy.concatenate(Overall_data['PixelError1'].to_numpy())
    print("Pooled performance",TrueFractionPixels.shape,PredfractionPixels.shape,IntensityPixels.shape,ErrorPixels.shape)
    
    fig, ax = plt.subplots(figsize=(8,6))
    #m3=ax.scatter(TrueFractionPixels, PredfractionPixels,c=IntensityPixels,cmap='turbo', label="STED Power = 30%")
    x_pred = numpy.linspace(0, 1, 10)
    lines=numpy.array([[Overall_data["fit_intercept1"][i],Overall_data["fit_slope1"][i],Overall_data["fit_correlation1"][i]] for i in range(len(pairs))])
    #print("LINES",lines.dtype)
    meanline=numpy.mean(lines.astype(float),axis=0)
    stdline=numpy.std(lines.astype(float),axis=0)
    negstdy=(meanline[0]-stdline[0])+ (meanline[1]-stdline[1]) * x_pred
    posstdy=(meanline[0]+stdline[0])+ (meanline[1]+stdline[1]) * x_pred
    meany=meanline[0]+meanline[1]*x_pred
    #for i in range(lines.shape[0]):
    #    line=lines[i,:]
    #    y_pred = line[0]+ line[1] * x_pred
        #ax.plot(x_pred, y_pred, label="Correlation r={}".format(line[2]))
    ax.fill_between(x_pred,negstdy,posstdy,alpha=0.3)
    ax.plot(x_pred,meany,label="Correlation r={:5.3f}+-{:5.3f}".format(meanline[2],stdline[2]))
    ax.plot([0,1],[0,1],'k--',label="Optimal")
    ax.set_xlabel('Ground Truth Fraction')
    ax.set_ylabel('Predicted Fraction')
    ax.set_ylim([0,1])
    ax.legend()
    
    counts, xbins, ybins = numpy.histogram2d(TrueFractionPixels,IntensityPixels, bins=(5,30))
    print("counts",counts.shape)
    print("xbins, ybins",xbins.shape, ybins.shape)
    sums, l, ll = numpy.histogram2d(TrueFractionPixels,IntensityPixels,weights=ErrorPixels, bins=(xbins, ybins))
    
    
    
    fig2,ax2=plt.subplots(figsize=(12,6))
    m3 =ax2.pcolormesh(ybins, xbins, sums / counts, cmap='Reds',vmin=0,vmax=1)
    plt.colorbar(m3, ax=ax2,label='Mean Error')
    ax2.set_xlabel('Intensity of pixel')
    ax2.set_ylabel('Fraction ')
    #ax2.legend()
    fig.savefig(os.path.join(savefolder,"Confusionlines_2species_fraction1.pdf"),transparent='True', bbox_inches="tight")
    fig2.savefig(os.path.join(savefolder,"ErrorColormap_2species_fraction1.pdf"),transparent='True', bbox_inches="tight")
    plt.close(fig)
    plt.close(fig2)

    TrueFractionPixels=numpy.concatenate(Overall_data["TrueFraction2"].to_numpy())
    PredfractionPixels = numpy.concatenate(Overall_data["PredictedFraction2"].to_numpy())
    IntensityPixels= numpy.concatenate(Overall_data["PixelIntensity2"].to_numpy())
    ErrorPixels=numpy.concatenate(Overall_data['PixelError2'].to_numpy())
    print("Pooled performance",TrueFractionPixels.shape,PredfractionPixels.shape,IntensityPixels.shape,ErrorPixels.shape)
    
    fig, ax = plt.subplots(figsize=(8,6))
    #m3=ax.scatter(TrueFractionPixels, PredfractionPixels,c=IntensityPixels,cmap='turbo', label="STED Power = 30%")
    x_pred = numpy.linspace(0, 1, 10)
    lines=numpy.array([[Overall_data["fit_intercept2"][i],Overall_data["fit_slope2"][i],Overall_data["fit_correlation2"][i]] for i in range(len(pairs))])
    #print("LINES",lines.dtype)
    meanline=numpy.mean(lines.astype(float),axis=0)
    stdline=numpy.std(lines.astype(float),axis=0)
    negstdy=(meanline[0]-stdline[0])+ (meanline[1]-stdline[1]) * x_pred
    posstdy=(meanline[0]+stdline[0])+ (meanline[1]+stdline[1]) * x_pred
    meany=meanline[0]+meanline[1]*x_pred
    #for i in range(lines.shape[0]):
    #    line=lines[i,:]
    #    y_pred = line[0]+ line[1] * x_pred
        #ax.plot(x_pred, y_pred, label="Correlation r={}".format(line[2]))
    ax.fill_between(x_pred,negstdy,posstdy,alpha=0.3)
    ax.plot(x_pred,meany,label="Correlation r={:5.3f}+-{:5.3f}".format(meanline[2],stdline[2]))
    ax.plot([0,1],[0,1],'k--',label="Optimal")
    ax.set_xlabel('Ground Truth Fraction')
    ax.set_ylabel('Predicted Fraction')
    ax.set_ylim([0,1])
    ax.legend()
    
    counts, xbins, ybins = numpy.histogram2d(TrueFractionPixels,IntensityPixels, bins=(5, 30))
    print("counts",counts.shape)
    print("xbins, ybins",xbins.shape, ybins.shape)
    sums, l, ll = numpy.histogram2d(TrueFractionPixels,IntensityPixels,weights=ErrorPixels, bins=(xbins, ybins))
    
    
    
    fig2,ax2=plt.subplots(figsize=(12,6))
    m3 =ax2.pcolormesh(ybins, xbins, sums / counts, cmap='Reds',vmin=0,vmax=1)
    plt.colorbar(m3, ax=ax2,label='Mean Error')
    ax2.set_xlabel('Intensity of pixel')
    ax2.set_ylabel('Fraction ')
    #ax2.legend()
    fig.savefig(os.path.join(savefolder,"Confusionlines_2species_fraction2.pdf"),transparent='True', bbox_inches="tight")
    fig2.savefig(os.path.join(savefolder,"ErrorColormap_2species_fraction2.pdf"),transparent='True', bbox_inches="tight")
    plt.close(fig)
    plt.close(fig2)
    Overall_data.to_csv(os.path.join(savefolder,"Overall_data_{}.csv".format(STEDPOWER)))
    

    return numpy.array([Overall_data["resolution1"],Overall_data["resolution2"],Overall_data["res_fraction1"],Overall_data["res_fraction2"],Overall_data["squirrel_f1"],Overall_data["squirrel_f2"]])
#["_30",[0,4]],["30",[0,0]],["40",[0,0]]

#[[5,[0,0]],[10,[0,0]],[15,[0,0]],[20,[0,0]],[30,[0,0]],[40,[0,0]]]
Powerslist=[[20,[0,0]]]
Powerslist=[[10,[0,0]],[20,[0,0]],[30,[0,0]],[40,[0,0]]]
#Powerslist=[[5,[0,0]],[10,[0,0]],[15,[0,0]],[20,[0,0]]]
globalcumstats=[]
globalcumstatsmean=[]
globalcumstatsstd=[]
row=0

for power in Powerslist:
    stats=Simulate2SpeciesSTED(power[0],power[1])
    print(stats.shape)

    plt.show()