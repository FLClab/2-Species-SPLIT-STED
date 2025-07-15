"""

Program that creates synthetic 2 species STED-FLIM images by summing single species images and then 
unmixing them using the 2 Species STED-FLIM method.
The program then compares the unmixed images to the ground truth images and computes different metrics such as resolution and nanoJ-SQUIRREL


The routine is defined as a function and called at the end of the script in a loop over the different STED powers.

"""
import os

from sys import path as path1;
Functionspath=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Functions")
path1.append(Functionspath)
from Main_functions import (line_equation, to_polar_coord, polar_to_cart, load_image,select_channel, get_foreground)
from Phasor_functions import Median_Phasor,DTCWT_Phasor,unmix2species
from objectives import (Squirrel, Bleach)
import decorr

from tiffwrapper import imsave,LifetimeOverlayer
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy
import glob
import itertools
import seaborn
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression
import skimage
from skimage import filters
import scipy

import easygui
import tifffile

#plt.style.use('dark_background')
matplotlib.rcParams['axes.linewidth'] = 0.8
# ------------------ Default Input variables ----------------
params_dict = {
    # Parameter in option in the matlab code
    "smooth_factor": 0.2,  # % 'Smoothing factor:'
    "phasor_smooth_cycles": 1,  # % 'Smoothing cycles phasor:'
    "foreground_threshold": 10,
    "harm1": 1,
}
# Define the folders containing the control images
filename1=easygui.diropenbox(default=os.path.expanduser("~Desktop"),title="Select folder containing control images for First fluorophore")
filename2=easygui.diropenbox(default=os.path.expanduser("~Desktop"),title="Select folder containing control images for Second fluorophore")

# List of powers to loop through and ID of the images to use as controls
Powerslist=[[10,[0,0]],[20,[0,0]],[30,[0,0]],[40,[0,0]]]
#Powerslist=[[5,[0,0]],[10,[0,0]],[15,[0,0]],[20,[0,0]]]

labels=["Bassoon_CF594",'PSD95_STORANGE',"Mixture"]


# Keys to the channels to use. For Tiff files, use channel ID
keyscontrols = ['STED 561 {11}','STED 561 {11}']
keysmixed = ['STED 561 {11}','STED 561 {11}']


#keyscontrols = ['STED_635P {2}', 'STED_635P {2}']
#keysmixed = [ 'STED_635P {2}', 'STED_635P {2}']
#keyscontrols=[1,1] # For Tiff file, use channel ID
#keysmixed=[1,1] # For Tiff file, use channel ID
# -----------------------------------------------------------
def Simulate2SpeciesSTED(STEDPOWER,NUMIM):
    """
    Generates simulated combinations of pairs of images acquired with the same STED power
      and performs 2species STED-FLIM to unmix them using control images designated by numim.

    Args:
        STEDPOWER (int): The STED power.
        NUMIM (list): A list of integers representing the control images to use for each species.

    Returns:
        An array containing the resolutions and nanoJ-SQUIRREL scores for all the pairs of images of the STED power.

    Saves : 
    For each pair of images:
        - A pdf file of the phasor plot colorcoded by the mixture fraction
        - ground truth and predicted fraction images
    Overall:
        - A CSV file of the data.
        - A legend file containing the paths to the control images and to the files that are mixed together for each pairID
    """

    f2= os.path.join(filename2,"*_{}percentSTED.msr".format(STEDPOWER))
    f1= os.path.join(filename1,"*_{}percentSTED.msr".format(STEDPOWER))
 
    filenames = [f1,f2]
    filenamescontrol = [f1,f2]
    
    
    colors=['royalblue','orangered','springgreen']
  
    # Create a folder to save the results
    savefolder = "Simulation_Cy3_{}Percent_2Species_PSD95Bassoon".format(STEDPOWER)
    savefolder = os.path.join(os.path.expanduser("~/Desktop"), savefolder)
    os.makedirs(savefolder,exist_ok=True)

    msrfiles = []   
    for k, filename in enumerate(filenamescontrol):
        print(labels[k])
        # Make list of all the images in the folder
        images = glob.glob(filename)
        # Select the image to use as control
        numim = NUMIM[k]
        image = images[numim]
        msrfiles.append(image)
    print(msrfiles)
    
    CoM_x, CoM_y = [], []
    # Create a figure for the phasor plot
    fig4,ax_scatter = plt.subplots(figsize=(3,3))
    # draw universal semi-circle
    edge = numpy.linspace(start=0, stop=15, num=200)
    theta = numpy.linspace(0, numpy.pi, 100)
    r = 0.5
    x1 = r * numpy.cos(theta) + 0.5
    x2 = r * numpy.sin(theta)
    ax_scatter.plot(x1, x2, color="k", ls="--",linewidth=0.8)
    
    # Create a legend file
    with open(os.path.join(savefolder,'legend.txt'),'w') as data: 
        data.write("Controls\n")
    scatterlist=[]

    for i, msr in enumerate(msrfiles) : 
        df = pd.DataFrame(columns=['x','y'])
        dg = pd.DataFrame(columns=['g', 's'])
        # Write the control images info in the legend file
        with open(os.path.join(savefolder,'legend.txt'),'a') as data: 
            data.write("{}\t{}\t{}\n".format(labels[i],keyscontrols[i],msr))

        # Load the image and select the channel    
        imagemsr = load_image(msr)
        image1= select_channel(imagemsr, keyscontrols[i])
        
        imsum = image1[:,:,10:111].sum(axis=2)
        imsum = imsum.astype('int16')

        # Calculate the phasor distribution of the foreground of the control image
        seuil = get_foreground(imsum)
        print("Caclulation for an image of shape", image1.shape, "...")
        params_dict["foreground_threshold"] = seuil
        
        print("foreground_threshold=", params_dict["foreground_threshold"])
    
        x, y, g_smoothed, s_smoothed, original_idxes = Median_Phasor(image1, params_dict, **params_dict)
        df['x'] = x.flatten()
        df['y'] = y.flatten()
        m, phi = to_polar_coord(df['x'], df['y'])
        g, s = polar_to_cart(m, phi)
        dg['g'], dg['s'] = g, s

        # Find the centroid of the phasor distribution using KMeans clustering
        kmeans = KMeans(n_clusters=1, init='k-means++', random_state=42)
        y_kmeans = kmeans.fit_predict(dg)
        CoM_x.extend(kmeans.cluster_centers_[:, 0][:].tolist())
        CoM_y.extend(kmeans.cluster_centers_[:, 1][:].tolist())

        # Plot the phasor distribution
        a=ax_scatter.scatter(g, s, s=1, c=colors[i], alpha=0.10,rasterized=True)
        scatterlist.append(a)


    ax_scatter.set_xlim(0, 1)
    ax_scatter.set_ylim(0, 1)

    ax_scatter.set_xlabel('g')
    ax_scatter.set_ylabel('s')
   
    xaxis = numpy.linspace(0, 1, 100)

    
    ##Calculate the endpoints of the linear trajectory using the centroids
    Pn_x, Pn_y = CoM_x[0], CoM_y[0]
    P_n = numpy.array([Pn_x, Pn_y])
    P2_x, P2_y = CoM_x[1], CoM_y[1]
    p2 = numpy.array([P2_x, P2_y])
    ## Line equation for the trajectory
    x2, y2 = [Pn_x, P2_x], [Pn_y, P2_y]
    m2, c2 = line_equation(x2, y2)
    y2 = m2 * xaxis + c2
    

    # Plot the trajectory line and the endpoints
    pnscatter = ax_scatter.scatter(Pn_x, Pn_y, s=50, c='gold')
    p2scatter = ax_scatter.scatter(P2_x, P2_y, s=50, c='gold')
    p2pnline = ax_scatter.plot([Pn_x, P2_x], [Pn_y, P2_y], c='dodgerblue')

    # Save the phasor plot
    fig4.savefig(os.path.join(savefolder, "Phasor_2species_ControlsOnly.pdf"), transparent='True', bbox_inches="tight")

    # Remove the previous scatter plots from the figure
    pnscatter.remove()
    p2scatter.remove()
    t = [scatter.remove() for scatter in scatterlist]
    ax_scatter.lines[-1].remove()



    # Combine all possible pairs of images, perform unmixing and calculate metrics
    msrfiles=[]
    images=[glob.glob(filename)for filename in filenames]
    number = [len(glob.glob(filename)) for filename in filenames]
    print('There are ',number, 'Images in these folders')
    pairs = list(itertools.product(images[0], images[1]))
    print(len(pairs))
    # Create a DataFrame to store the results
    Overall_data = pd.DataFrame(
        columns=["Power",'image1', 'image2', 'resolution1', 'resolution2',"res_fraction1", "res_fraction2","squirrel_f1","squirrel_f2",
                 "squirrelsmooth_f1","squirrelsmooth_f2"])
    # Loop through each pair of images
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
        # For each image in the pair, calculate the resolution and the foreground mask
        for i, msr in enumerate(msrfiles):
            print("i",i)
            imagemsr=load_image(msr)
       
            image1 = select_channel(imagemsr, keysmixed[i])
            imagec1 = select_channel(imagemsr, keyscontrols[i])

            res_mix = decorr.calculate(numpy.sum(image1[:, :, 10:111],axis=2))
            if math.isinf(res_mix):
                res_mix=10
      
            ov_data.append(res_mix * 20)
       
            Imagelist.append(image1)
            print(image1.shape)
            imsum=numpy.sum(image1[:,:,10:111],axis=2)

            seuil=3
            seuils.append(seuil)
            mask=imsum>seuil
            mask=scipy.ndimage.binary_fill_holes(mask)
            Masklist.append(mask)
    
        print("masklist",len(Masklist),Masklist[0].shape)
         # Create empty images of the correct size to store the combined image and masks
        Combo=numpy.zeros(numpy.max([Imagelist[0].shape,Imagelist[1].shape],axis=0))
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
        # Add the first image to the combined image 
        Combo[minx:maxx,miny:maxy ,:]+=Imagelist[0][:,:, :]
        Combomask1[minx:maxx,miny:maxy]+=Masklist[0]
    
        Combosingle1[minx:maxx, miny:maxy] += numpy.sum(Imagelist[0][:,:,10:111],axis=2)
        minx = 0
        maxx = Imagelist[1].shape[0]
        miny = 0
        maxy = Imagelist[1].shape[1]
        print(minx, miny, maxx, maxy)
        # Add the second image to the combined image 
        Combo[minx:maxx, miny:maxy, :] += Imagelist[1][:, :, :]
        Combomask2[minx:maxx, miny:maxy] += Masklist[1]
    
        Combosingle2[minx:maxx, miny:maxy] += numpy.sum(Imagelist[1][:,:,10:111], axis=2)
        ComboMasklist = [Combomask1, Combomask2]
        ComboSinglelist=[Combosingle1,Combosingle2]
    
        Imagelist.append(Combo)
        
        image1=Combo.copy()
        # Calculate the phasor distribution of the foreground of the combined image
        
        CoM_x, CoM_y = [], []
        d_melange = pd.DataFrame(columns=['g', 's'])
        df = pd.DataFrame(columns=['x', 'y'])
        dg = pd.DataFrame(columns=['g', 's'])
        print(image1.shape)
        imsum = image1[:,:,10:111].sum(axis=2)
        imsum = imsum.astype('int16')
        seuil=min(seuils)
        print("Caclulation for an image of shape", image1.shape, "...")
        params_dict["foreground_threshold"] = seuil
        print("foreground_threshold=", params_dict["foreground_threshold"])
        #x,y, original_idxes,Images,Images_Filtered=DTCWT_Phasor(image1, 0, nlevels=10, neighborhood=50)
        x, y, g_smoothed, s_smoothed, original_idxes = Median_Phasor(image1, params_dict, **params_dict)
        #x = x[imsum > params_dict["foreground_threshold"]]
        #y = y[imsum > params_dict["foreground_threshold"]]
        df['x'] = x.flatten()
        df['y'] = y.flatten()
        # Calibrate the phasor distribution using the IRF measurement
        m, phi = to_polar_coord(df['x'], df['y'])
        g, s = polar_to_cart(m, phi)
        dg['g'], dg['s'] = g, s

        p3 = dg[['g', 's']].to_numpy() #phasor of the combined image

        # Perform unmixing of the combined image using the 2 species STED-FLIM method
        imsum_flat_lin1, imsum_flat_lin2, Solve = unmix2species(p3, original_idxes, Imagelist[2], P_n, p2)
    
        fraction2 = imsum_flat_lin1.copy()
        fraction1 = imsum_flat_lin2.copy()
    
    # Calculate the fraction of each species in the mixture in the ground truth and in the prediction
        GroundTruth_Fraction=[(Combosingle1/(Combosingle2+Combosingle1)),(Combosingle2 / (Combosingle2 + Combosingle1))]
        Predicted_Fraction=[fraction2.copy(),fraction1.copy()]


# Plot the color-coded phasor plot and the line used for unmixing
        mixphasor = ax_scatter.scatter(g, s, s=2,c=Solve[1,:],cmap="cool",rasterized=True,label="Mixture")
        lineplot=ax_scatter.plot(xaxis, y2, 'dodgerblue')
        pnplot=ax_scatter.scatter(Pn_x, Pn_y, s=50,c='gold')
        p2plot=ax_scatter.scatter(P2_x, P2_y, s=50,c='gold')
        # Save the phasor plot
        fig4.savefig(os.path.join(savefolder,"Phasor_2species_{}.pdf".format(Pair_id)),transparent='True', bbox_inches="tight",dpi=900)
        # Remove the previous scatter plots from the figure
        pnplot.remove()
        p2plot.remove()
        ax_scatter.lines[-1].remove()

        mixphasor.remove()
        
    
        Combomask=numpy.logical_or(ComboMasklist[0],ComboMasklist[1])

         # Apply the calculated fractions to the intensity of the combined image to get the fraction images
        imsum = Combo[:, :, 10:111].sum(axis=2)
        fraction1 *= imsum
        fraction2 *= imsum
        # Calculate the resolution of the fractions
        res_fraction1 = decorr.calculate(fraction1.astype(numpy.uint16) )
        res_fraction2 = decorr.calculate(fraction2.astype(numpy.uint16) )
        # If the resolution does not converge, set it to a default value
        if math.isinf(res_fraction2):
            res_fraction2=10
        if math.isinf(res_fraction1):
            res_fraction1=10
        print("res_fraction1", res_fraction1 * 20, "res_fraction2", res_fraction2 * 20)


        # Calculate the NanoJ-SQUIRREL scores and error maps for the fractions
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

        # Add the resolutions and the NanoJ-SQUIRREL scores to the dataframe
        ov_data.extend([res_fraction1* 20, res_fraction2* 20,squirrel_f1,squirrel_f2,squirrelsmoothf1[2],squirrelsmoothf2[2]])
        Overall_data.loc[Pair_id] = ov_data

        # Save the SQUIRREL error maps to tiff files
        imagecomp=numpy.dstack((squirrelmap1,squirrelmap2))
        imagecomp=numpy.moveaxis(imagecomp,2,0)
        filenameout = os.path.join(savefolder,"{}_SquirrelMaps.tiff".format(Pair_id))
        tifffile.imwrite(filenameout, imagecomp)

        # Save the ground truth and predicted fraction images to a tiff file
        imagecomp=numpy.dstack((GroundTruth_Fraction[0],GroundTruth_Fraction[1],Combosingle1,Combosingle2,imsum))
        imagecomp=numpy.moveaxis(imagecomp,2,0)
        filenameout =os.path.join(savefolder, "{}_GroundTruth.tiff".format(Pair_id))
        print(filenameout)
        tifffile.imwrite(filenameout, imagecomp)

        # Save the individual masks and the combined mask to a tiff file
        imagecomp=numpy.dstack((ComboMasklist[0],ComboMasklist[1],Combomask))
        imagecomp=numpy.moveaxis(imagecomp,2,0)
        filenameout = os.path.join(savefolder,"{}_Masks.tiff".format(Pair_id))
        tifffile.imwrite(filenameout, imagecomp)

        # Save the unmixed fraction images to a tiff file
        imagecomp=numpy.dstack((Predicted_Fraction[0],Predicted_Fraction[1],fraction1,fraction2))
        imagecomp=numpy.moveaxis(imagecomp,2,0)
        filenameout = os.path.join(savefolder,"{}_Predictions.tiff".format(Pair_id))
        tifffile.imwrite(filenameout, imagecomp)

    # Create a composite image of the ground truth images and save it to a tiff file
        imagecomp=numpy.dstack((Combosingle1,Combosingle2))
        imagecomp=numpy.moveaxis(imagecomp,2,0)
        filenameout =os.path.join(savefolder, "{}_GroundTruth_Composite.tiff".format(Pair_id))
        imsave(file=filenameout, data=imagecomp.astype(numpy.uint16), composite=True, luts=("Cyan Hot","Magenta Hot"), pixelsize=(20E-3,20E-3))
      # Create a composite image of the unmixed images and save it to a tiff file
        imagecomp=numpy.dstack((fraction2,fraction1))
        imagecomp=numpy.moveaxis(imagecomp,2,0)
        filenameout =os.path.join(savefolder, "{}_Predicted_Composite.tiff".format(Pair_id))
        imsave(file=filenameout, data=imagecomp.astype(numpy.uint16), composite=True, luts=("Cyan Hot","Magenta Hot"), pixelsize=(20E-3,20E-3))
       
        with open(os.path.join(savefolder,'legend.txt'),'a') as data: 
            data.write("{}\t{}\t{}\n".format(Pair_id,a,b))
    

        # Make a figure to display the ground truth and predicted fraction images and the error maps
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
    

        # Create an overlay of the unmixed fraction image on the intensity images
        overlayer = LifetimeOverlayer(fraction1, imsum/imsum.max(), cname='CET-I1')
        lifetime_rgb, cmap = overlayer.get_overlay(
             lifetime_minmax=(0., 1),
             intensity_minmax=(0, 0.5)  # inTensity saturated to get more bright regions
         )
        filenameout = os.path.join(savefolder,"{}_OverlayF1.tiff".format(Pair_id))
        tifffile.imwrite(filenameout, lifetime_rgb.astype(numpy.float32))    

    
    Overall_data.to_csv(os.path.join(savefolder,"Overall_data_{}.csv".format(STEDPOWER)))
    

    return numpy.array([Overall_data["resolution1"],Overall_data["resolution2"],Overall_data["res_fraction1"],Overall_data["res_fraction2"],Overall_data["squirrel_f1"],Overall_data["squirrel_f2"]])


# ------------------ Main Simulation Loop ----------------

for power in Powerslist:
    stats=Simulate2SpeciesSTED(power[0],power[1])
    print(stats.shape)

    plt.show()