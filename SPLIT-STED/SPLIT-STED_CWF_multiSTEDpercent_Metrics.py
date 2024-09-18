
"""
Takes a folder of microscope files containing 4 images (Confocal (intensity only), Confocal-FLIM, STED-FLIM, Confocal (intensity only)))
Computes :
    SPLIT-STED image with DTCWT filtering from the STED-FLIM image
    Resolution of Confocal FLIM, STED FLIM ans SPLIT STED images
    SQUIRREL metrics of SPLIT-STED compared to STED
    Photobleaching caused by acquisition of STED-FLIM image
Outputs:
    - Dataframe with all the metrics
    - Plots of the resolution of the images
    - Plots of the phasor space of the images
    - Tiff image files of Confocal, STED and SPLIT-STED images
    - Fractional component maps of SPLIT-STED image

"""

from skimage import io
import os
import glob
import numpy

import tifffile
import easygui
import skimage.io as skio
import numpy as np
import scipy 
from scipy.optimize import curve_fit
from pandas.plotting import table
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


import os
import easygui
from plotly.offline import plot


from sklearn.cluster import KMeans

import os.path
from sys import path as syspath; 
dossier = os.path.expanduser("~/Documents/Github/Abberior-STED-FLIM/Functions")
syspath.append(dossier)
from decorr_res import decorr_res
from objectives import (Squirrel, Bleach)
from Main_functions import (line_equation, to_polar_coord, polar_to_cart, choose_msr_file, get_foreground)
from Phasor_functions import DTCWT_Phasor,Median_Phasor,SPLIT_STED
from matplotlib.gridspec import GridSpec
from convertmsr_bioformatsAB import MSRReader
from tiffwrapper import imsave,LifetimeOverlayer
import time
import matplotlib
import sklearn
matplotlib.rcParams['axes.linewidth'] = 0.8
# -----------------------------------------------------------
#    Sélection des images dans un même fichier avec easygui
#easygui.diropenbox(default=os.path.expanduser("~Desktop"))
params_dict = {

    # Parameter in option in the matlab code
    #    "Tg" : 6, #% 'First frame to sum:'
    "Nb_to_sum": 250,  # The Tg infered from this variable override Tg
    "smooth_factor": 0.08,  # % 'Smoothing factor:'
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
nlevels=2
neighbours=50

#plt.style.use('dark_background')

## Path of folder containing images

#filename = easygui.diropenbox(default=os.path.expanduser("~Desktop"))
filename = os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5',"PSD95_AF647_STEDPowerBleach_5to20_2")
#filename=os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5','rab_Bassoon_STAR635P_STEDPowerBleach_5to20_2')
filename= os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5',"msB2Spectrin_AF647_STEDPowerBleach_5to30_1")
#filename= os.path.join('U:', os.sep,'adeschenes','2024-02-29_FLIM_Cy5',"alphaTubulin_AF647_STEDPowerBleach_5to20_1")
#
filename= os.path.join('U:', os.sep,'adeschenes','2024-03-06_FLIM_PSDBassoon_Cy3',"AlphaTubulin_AF647_1to500_STEDPowerBleach_1")
filename= os.path.join('T:', os.sep,'adeschenes','SimulationDataset_STEDFLIM','Cy3',"PSD95_STORANGE")
#filename= os.path.join('T:', os.sep,'adeschenes','SimulationDataset_STEDFLIM','Cy3',"rabBassoon_CF594")
filename  = os.path.join('T:', os.sep, 'adeschenes', 'SimulationDataset_STEDFLIM', 'Cy3', 'Homer_STORANGE',"MediumAcq")
filename= os.path.join('T:', os.sep,'adeschenes','SimulationDataset_STEDFLIM','Cy5','alphaTubulin_Alexa647')
filename= os.path.join('T:', os.sep,'adeschenes','SimulationDataset_STEDFLIM','Cy5','B2Spectrin Alexa647')
filename= os.path.join('T:', os.sep,'adeschenes','SimulationDataset_STEDFLIM','Cy5','rabBassoon STAR635P')
## Dictionary of the images to be used (keys in the .msr files)

# mapcomp = { 'Conf pre'  :  'Conf_Pre {13}',
#               'Conf FLIM' :  'Conf640 {10}',
#             'STED FLIM' :'STED640 {10}',
#             'Conf post' : 'Conf_Post {14}',
#            'STED High' : 'STED_635P_HighP {8}'}

# mapcomp = { 'Conf pre'  : 'Confocal_Pre {14}',
#              'Conf FLIM' : 'Confocal_561 {11}',
#              'STED FLIM' :'STED 561 {11}',
#              'Conf post' : 'Confocal_Post {15}',
#             'STED High' : 'STED 561_HighP {16}'}

mapcomp = { 'Conf pre'  : 'Conf_pre {6}',
             'Conf FLIM' : 'Conf_635P {2}', 
             'STED FLIM' :'STED_635P {2}',
             'Conf post' :  'Conf_post {7}',
            'STED High' : 'STED 561_HighP {16}'}


colors=["springgreen",'orangered','gold','deepskyblue']
labels=["STED Phasor","Phasor centroids","Pure species coordinates","Projection Line"]

# Path of the folder where the results will be saved
savefolder=str(input("Name of Output folder: "))
savefolder = os.path.join(os.path.expanduser("~/Desktop"), savefolder)
os.makedirs(savefolder, exist_ok=True)

# List of images in the folder
path=os.path.join(filename,"*.msr")
images = glob.glob(path)
print("path",path)
numim=len(images)
print('There are ',numim, 'Images in this folder')

# Create a dataframe to store the metric measurements
im_ids = []
data_objectifs = pd.DataFrame(columns=['Filename',
                'STEDpercent','Res conf stack', 'Res sted stack', 'Res splitsted', 'Squirrel','Bleach','res_HighSTED',"Squirrel HighP","Squirrel SPLITvsConf"])

## Create colormaps for color-coded phasors. Magenta to yellow  and cyan to yellow 
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    name="coolspring",
    colors=["#00ffffff", "#ffff00ff"]
)

matplotlib.colormaps.register(cmap=cmap, force=True)
matplotlib.colormaps.register(cmap=cmap.reversed(), force=True)

cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    name="springcool",
    colors=["#ff00ffff", "#ffff00ff"]
)
matplotlib.colormaps.register(cmap=cmap, force=True)
matplotlib.colormaps.register(cmap=cmap.reversed(), force=True)

phasorcolor="springcool"
#phasorcolor="springcool"
with MSRReader() as msrreader:
# Create a figure to plot the phasors
    fig4,ax4 = plt.subplots(figsize=(2,2))
    ax4.set_xlim(0, 1)
    ax4.set_ylim(0, 1)

    fig5,ax5 = plt.subplots(figsize=(2,2))
    ax5.set_xlim(0, 1)
    ax5.set_ylim(0, 1)

    ax4.set_xlabel('g')
    ax4.set_ylabel('s')
# Create universal semi-circle for phasor space and add it to the plot
    edge = np.linspace(start=0, stop=15, num=200)
    theta = np.linspace(0, np.pi, 100)
    r = 0.5
    x1 = r * np.cos(theta) + 0.5
    x2 = r * np.sin(theta)
    ax4.plot(x1, x2, color="black", ls="--",linewidth=0.8)
    ax5.plot(x1, x2, color="black", ls="--", linewidth=0.8)

    # Loop through images in the folder
    for i,im in enumerate(images) : 
        print("######################")
        print(i,"of",numim, os.path.basename(im))
        print("######################")
        imagemsr = msrreader.read(im)
        print(imagemsr.keys())
        print('image opened with success')

        #Extract the depletion power from the image filename
        image_id=i
        sted_percent = str(os.path.basename(im).split('_')[-1].split('percentSTED')[0])

        #Extract the intensity confocal images from the image file
        Conf_init = imagemsr[mapcomp['Conf pre']]
        Conf_end = imagemsr[mapcomp['Conf post']]
        objec = []
        objec.append(os.path.basename(im))
        objec.append(sted_percent)
        fg = []

        for i,key in enumerate(mapcomp) : 

            if key == 'Conf pre' :
                continue
            if key =='Conf post' :
                continue
            if key == 'STED High' :
                if mapcomp['STED High'] in imagemsr.keys():   
                    image1 = imagemsr[mapcomp[key]]
                    imsumhigh = image1[:,:,10:111].sum(axis=2)

                    filenameout = os.path.join(savefolder,
                                           os.path.basename(im).split(".msr")[0] + "_HighP_STED.tiff")
                    imsave(filenameout, imsumhigh.astype(numpy.uint16), luts="Red Hot")
                    res_HighSTED = decorr_res(imname=None, image=imsumhigh)
                    objec.append(res_HighSTED)
                    fg_highpSTED=get_foreground(imsumhigh)   
                    squirrel_highp_vs_conf = Squirrel(method="L-BFGS-B", normalize=True).evaluate([imsumhigh], conf_stack, conf_stack,imsumhigh > fg_highpSTED, conf_stack >fg_conf_stack )
                    squirrel_sted_vs_conf = Squirrel(method="L-BFGS-B", normalize=True).evaluate([im_splitsted], sted_stack, conf_stack,im_splitsted > fg_splitsted, sted_stack >fg_sted_stack )
                    objec.extend([squirrel_highp_vs_conf,squirrel_sted_vs_conf])
                else:
                    print("No high power images taken")
                    objec.extend([numpy.nan,numpy.nan,numpy.nan])
                continue
            CoM_x, CoM_y = [], []
            
            df = pd.DataFrame(columns=['x','y'])
            dg = pd.DataFrame(columns=['g', 's'])
            df1 = pd.DataFrame(columns=['x','y'])
            dg1 = pd.DataFrame(columns=['g', 's'])

# Measure the foreground threshold of the image
            image1 = imagemsr[mapcomp[key]]
            imsum = image1[:,:,10:111].sum(axis=2)

            seuil = get_foreground(image1)
            fg.append(seuil)
            #print("Caclulation for an image of shape", image1.shape, "...")
            params_dict["foreground_threshold"] = seuil

            params_dict["Nb_to_sum"] = image1.shape[2]
            print("foreground_threshold=", params_dict["foreground_threshold"])


# Calculate the phasor distribution with DTCWT filtering
            start_DTCWT_Time = time.time()
            x, y, original_indexs, Images, Images_Filtered = DTCWT_Phasor(image1, 0, nlevels,neighbours)
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
            g,s = numpy.array(g),numpy.array(s)
            indexes = numpy.where((g > 0) & (g < 1) & (s > 0) & (s < 1))
            g,s = g[indexes],s[indexes]
            original_idxes = Image_indices[indexes]
            #original_idxes = Image_indices

            dg['g'] = g
            dg['s'] = s

            
            if key == 'Conf FLIM' :
            # Save the confocal-flim image (sum intensity over time bins) and measure its resolution
                conf_stack = imsum
                filenameout = os.path.join(savefolder,
                                           os.path.basename(im).split(".msr")[0] + "_Confocal.tiff")
                imsave(filenameout, imsum.astype(numpy.uint16), luts="Red Hot")
                #tifffile.imwrite(filenameout, imsum.astype(numpy.uint16))
                res_conf = decorr_res(imname=None, image=imsum)
                objec.append(res_conf)
            # Find the centroid coordinates of the phasor distribution
                n=1
                kmeans = KMeans(n_clusters = n, init = 'k-means++', random_state = 42)
                y_kmeans = kmeans.fit_predict(dg)
                CoM_x, CoM_y = kmeans.cluster_centers_[:, 0], kmeans.cluster_centers_[:, 1]
            # Project the centroids on the semi-circle in the phasor space
                r = 0.5 
                norm = numpy.sqrt((CoM_x[0] - 0.5) ** 2 + (CoM_y[0] ** 2) )
                Pn_x = 0.5 + (r * (CoM_x[0] - 0.5) / norm)
                Pn_y = 0 + r * (CoM_y[0] - 0) / norm

                xaxis = np.linspace(0, 1, 100)

            #Line between P_n and (1,0)
                x2, y2 = numpy.asarray([Pn_x, 1.0], dtype = 'float64'), numpy.asarray([Pn_y, 0.0], dtype = 'float64')
                m2, c2 = line_equation(x2, y2)
                y2 = m2 * xaxis + c2
                P_n = numpy.array([Pn_x, Pn_y])
                p2 = numpy.array([1,0])

            # Plot the Confocal-FLIM phasor distribution 
                confphasor=ax5.scatter(g,s,c="lime",s=0.5,alpha=0.2,rasterized=True)


                continue
            

            #print(key, 'Im here now')
            sted_stack = imsum

        # Measure the resolution of the STED-FLIM image (intensity sum over time bins)
            res_sted_stack = decorr_res(imname=None, image=image1[:,:,10:].sum(axis=2))
            if numpy.isinf(res_sted_stack):
                res_sted_stack = 0
            objec.append(res_sted_stack)

  
        # Calculate the position of 2 centroids for STED-FLIM phasor distribution
            kmeans = KMeans(n_clusters = 2, init = 'k-means++', random_state = 42)
            y_kmeans = kmeans.fit_predict(dg)

            CoM_x.extend(kmeans.cluster_centers_[:, 0][:].tolist())
            CoM_y.extend(kmeans.cluster_centers_[:, 1][:].tolist())

            CoM_x.sort()
            CoM_y.sort(reverse=True) # To get the order CoM = [STED1 < STED2, ...]
        # Calculate the phasor of the STED-FLIM image with a lower foreground threshold to include all pixels in SPLIT-STED 
            start_DTCWT_Time = time.time()
            x, y, original_indexs, Images, Images_Filtered = DTCWT_Phasor(image1, 0,nlevels,neighbours)
            stop_DTCWT_Time = time.time()
            DTCWT_Time = stop_DTCWT_Time - start_DTCWT_Time
            params_dict["foreground_threshold"] = 3
            print("DTCWT_Time", DTCWT_Time)
            x = x[imsum > params_dict["foreground_threshold"]]
            y = y[imsum > params_dict["foreground_threshold"]]
            Image_indices = np.arange(len(imsum.flatten())).reshape(imsum.shape)
            Image_indices = Image_indices[imsum > params_dict["foreground_threshold"]]
        # Apply phasor calibration in polar coordinates (based on IRF measurement) and return to cartesan (g,s)
            df1['x'],df1['y'] = x.flatten(),y.flatten()
            x,y= x.flatten(),y.flatten()
            m, phi = to_polar_coord(x, y)
            g, s = polar_to_cart(m, phi)
            g,s = numpy.array(g),numpy.array(s)
            indexes = numpy.where((g > 0) & (g < 1) & (s > 0) & (s < 1))
            g,s = g[indexes], s[indexes]
            original_idxes = Image_indices[indexes]
            dg1['g'],dg1['s'] = g,s

            p3 = dg1[['g', 's']].to_numpy()

            p3_min = numpy.array([CoM_x[0], CoM_y[0]]) #STED1
            p3_max =numpy.array([CoM_x[1], CoM_y[1]]) #STED2
        # Plot the phasor distribution of the STED-FLIM image
            stedphasor=ax5.scatter(g, s, c="orange", s=0.5,alpha=0.2, rasterized=True)


            
        # Calculate the SPLIT-STED fractional components and apply to the STED-FLIM image
            im_fract, projection,t2 = SPLIT_STED(p3, P_n, p2, p3_min, p3_max, image1, original_idxes)

            im_splitsted = im_fract * imsum

            l2 = numpy.sum((P_n - p2) ** 2)  # distance between P_n and p2
            t_min = numpy.sum((p3_min - P_n) * (p2 - P_n)) / l2  #Project phasor centroids on line connecting p2 and P_n
            t_max = numpy.sum((p3_max - P_n) * (p2 - P_n)) / l2
            projectionp3min = P_n + numpy.multiply(p2 - P_n, t_min) # Find g,s coordinates of projected points (for phasor space graphs)
            projectionp3max = P_n + numpy.multiply(p2 - P_n, t_max) # Find g,s coordinates of projected points (for phasor space graphs)


        # Plot the color-coded phasor distribution of the SPLIT-STED image and line connecting the phasor centroids
            mixphasor = ax4.scatter(g, s, c=t2[1,:],cmap=phasorcolor, s=2,rasterized=True)
            p3minscat=ax4.scatter(projectionp3min[0],projectionp3min[1],s=50,c='orangered')
            p3maxscat=ax4.scatter(projectionp3max[0],projectionp3max[1],s=50,c='orangered')
            p2pnscat=ax4.scatter([P_n[0],p2[0]],[P_n[1],p2[1]],s=50,c='gold')
            lineplot=ax4.plot([P_n[0],p2[0]],[P_n[1],p2[1]],linewidth=3,c='deepskyblue')

        # Fit ellipse on the 70th percentile of the phasor distribution of the STED-FLIM image and plot it
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
                ax4.add_patch(ell)
                ell.set_facecolor("None")
                #ell.set_edgecolor(colors[k][a])
                ell.set_edgecolor("k")
                ell.set_linewidth(1.0)

# Save the SPLIT-STED Phasor graph
            fig4.savefig(os.path.join(savefolder,"Phasor_SPLITSTED_DTCWT_{}.pdf".format(os.path.basename(im).split(".msr")[0])),transparent='True', bbox_inches="tight",dpi=900)
 # Save the Confocal-FLIM and STED-FLIM Phasor graph and clean the plot for the next image       
            fig5.savefig(
                os.path.join(savefolder, "Phasor_raw_DTCWT_{}.pdf".format(os.path.basename(im).split(".msr")[0])),
                transparent='True', bbox_inches="tight", dpi=900)
            stedphasor.remove()
            confphasor.remove()
            p3minscat.remove()
            p3maxscat.remove()
            p2pnscat.remove()
            ax4.lines.pop(-1)
            mixphasor.remove()
            ell.remove()

    # Save the STED image, SPLIT-STED image and its fractional component maps     
            overlayer = LifetimeOverlayer(1-im_fract, imsum/imsum.max(), cname=phasorcolor)
            lifetime_rgb, cmap = overlayer.get_overlay(
            lifetime_minmax=(0, 1.0),
            intensity_minmax=(0, 0.6) # inensity saturated to get more bright regions
            )
            fig, ax = plt.subplots()
            ax.axis('off')
            ax.imshow(lifetime_rgb)
            cbar = fig.colorbar(cmap, ax=ax)
            fig.savefig(os.path.join(savefolder,os.path.basename(im).split(".msr")[0] + "_F1Overlay.pdf"),transparent='True', bbox_inches="tight")
            plt.close(fig)

            im_fractnan=im_fract.copy()
            im_fractnan[imsum < params_dict["foreground_threshold"]]=numpy.nan

            fig, ax = plt.subplots()
            ax.axis('off')
            ims=ax.imshow(im_fractnan,cmap=phasorcolor)
            fig.colorbar(ims, ax=ax)
            filenameout = os.path.join(savefolder,os.path.basename(im).split(".msr")[0] + "_SPLIT_STED_fractmap_{}_{}.pdf".format(nlevels,neighbours))
            #imsave(filenameout, im_fractnan, luts=phasorcolor)
            fig.savefig(filenameout,transparent='True', bbox_inches="tight")
            plt.close(fig)

            reversefract=(1-im_fract)
            reversefract[imsum < params_dict["foreground_threshold"]]=numpy.nan

            fig, ax = plt.subplots()
            ax.axis('off')
            ims=ax.imshow(reversefract,cmap=phasorcolor)
            fig.colorbar(ims,ax=ax)

            filenameout = os.path.join(savefolder,os.path.basename(im).split(".msr")[0] + "_SPLIT_STED_fractmap_reversed_{}_{}.pdf".format(nlevels,neighbours))
            #imsave(filenameout, reversefract, luts=phasorcolor)
            fig.savefig(filenameout,transparent='True', bbox_inches="tight")
            plt.close(fig)


            filenameout = os.path.join(savefolder,os.path.basename(im).split(".msr")[0] + "_SPLIT_STED_fractmap_{}_{}.tiff".format(nlevels,neighbours))
            imsave(filenameout, im_fract, luts="Red Hot")

            filenameout = os.path.join(savefolder,os.path.basename(im).split(".msr")[0] + "_SPLIT_STED_dtcwt_{}_{}.tiff".format(nlevels,neighbours))
            imsave(filenameout, im_splitsted.astype(numpy.uint16), luts="Red Hot")
            #tifffile.imwrite(filenameout, im_splitsted.astype(numpy.uint16))
            filenameout = os.path.join(savefolder,os.path.basename(im).split(".msr")[0] + "_STED.tiff")
            imsave(filenameout,image1[:,:,10:111].sum(axis=2).astype(numpy.uint16), luts="Red Hot")
            #tifffile.imwrite(filenameout, image1[:,:,10:].sum(axis=2).astype(numpy.uint16))
           
# Measure the resolution of the SPLIT-STED image
            res_splitsted = decorr_res(imname=None, image=im_splitsted)
            if numpy.isinf(res_splitsted):
                res_splitsted = 0
            objec.append(res_splitsted)

            fg_splitsted = get_foreground(im_splitsted)
            fg_sted_stack=get_foreground(sted_stack)
            fg_conf_stack=get_foreground(conf_stack)
   # Calculate the Squirrel metrics of the SPLIT-STED image compared to the STED-FLIM image         
            y_result = Squirrel(method="L-BFGS-B", normalize=True).evaluate([sted_stack], conf_stack, conf_stack, im_splitsted > fg_splitsted, conf_stack > fg_conf_stack )
            objec.append(y_result)
 # Calculate the photobleaching caused by the acquisition of the STED-FLIM image
            bleach = Bleach().evaluate(sted_stack , Conf_init, Conf_end, sted_stack > fg_sted_stack, Conf_init > fg_conf_stack)
   
            objec.append(bleach)
        if 0 not in objec[3:5]: #Don't include image in stats if resolution did not converge
            data_objectifs.loc[image_id] = objec
        else:
            print("Resolution didn't converge")
            print(res_splitsted ,res_sted_stack)
        plt.close(fig)
        
        

# Save the dataframe with the metrics and calculate the mean and standard deviation of the metrics
data_objectifs.STEDpercent= data_objectifs.STEDpercent.astype(float)
data_objectifs=data_objectifs.sort_values(by=['STEDpercent'])
data_objectifs.to_csv(os.path.join(savefolder,"data_objectifs.csv"))
print(data_objectifs)
dfmeanoverall=data_objectifs.mean()
dfstdoverall=data_objectifs.std()
dfmean=data_objectifs.groupby('STEDpercent', as_index=False).mean()
dfmean=dfmean.sort_values(by=['STEDpercent'])
print(dfmean[['Bleach']])
dfstd=data_objectifs.groupby('STEDpercent', as_index=False).std()
#print(dfstd)
dfstd=dfstd.sort_values(by=['STEDpercent'])
print(dfstd[['Bleach']])

# Plot the resolution of the images

Fig,Ax=plt.subplots(figsize=(12,8))
#Ax2 = Ax.twinx()
Ax.errorbar(x=dfmeanoverall['STEDpercent']*0,y=dfmeanoverall['Res conf stack']*20,yerr=dfstdoverall['Res conf stack']*20, fmt="o",c='lightskyblue',label='Confocal',ecolor='lightskyblue',capsize=10,elinewidth=3,ms=25)
#Ax.axhline(y=dfmeanoverall['Res conf stack']*20,ls='--',label='Confocal',color='lightskyblue',lw=3)
Ax.errorbar(x=dfmean['STEDpercent'],y=dfmean['Res sted stack']*20,yerr=dfstd['Res sted stack']*20, fmt="o",c='deepskyblue',label='STED',ecolor='deepskyblue',capsize=10,elinewidth=2,ms=25)
Ax.errorbar(x=dfmean['STEDpercent'],y=dfmean['Res splitsted']*20,yerr=dfstd['Res splitsted']*20, fmt="o",c= 'mediumblue',label='SPLIT-STED',ecolor='mediumblue',capsize=10,elinewidth=2,ms=25)

x=[data_objectifs['STEDpercent'][i] for i in range(data_objectifs['STEDpercent'].shape[0])]
#print("x",x)
xconf=[0 for i in range(data_objectifs['STEDpercent'].shape[0])]
x.extend(xconf)

#print("x",x)
y=[data_objectifs['Res sted stack'][i]*20 for i in range(data_objectifs['STEDpercent'].shape[0])]
#print("y",y)
yconf=[data_objectifs['Res conf stack'][i]*20 for i in range(data_objectifs['STEDpercent'].shape[0])]
y.extend(yconf)
pp=[dfmean['STEDpercent'][i] for i in range(dfmean['STEDpercent'].shape[0])]
#print("xy",x,y)
#xx=dfmeanoverall['STEDpercent']*0
#yy=dfmeanoverall['Res conf stack']*20
#print("xxyy",xx,yy)
#x=numpy.stack((x,xx))
#y=numpy.stack((y,yy))
def STEDRes(x, a, b):
    """
    Function to fit to resolution vs STED power graphs.

    Parameters:
    x (float): The STED power.
    a (float): Parameter alpha.
    b (float): Parameter dc.

    Returns:
    float: The calculated resolution.

    """
    return a/np.sqrt(1+(a**2*b**2*x))
xs = np.linspace(np.min(x), np.max(x),(np.max(x)-np.min(x)).astype(int))
# Fit the resolution vs STED power graph
try:

    params, cv = curve_fit(STEDRes,x,y,p0=(dfmeanoverall['Res conf stack']*20,3.4E-3),maxfev = 10000)
    m,b = params
    print("x",x[0].dtype,x)
    print("m",m.dtype,m)
    print("b",b.dtype,b)
    y_pred =STEDRes(np.array(x), m,b)
    rsquared=sklearn.metrics.r2_score(y, y_pred)
    Ax.plot( xs, STEDRes(xs,m,b), '-',c='deepskyblue')
    Ax.text(5,200,"Parameters: \nalpha={0:.3g} \ndc={1:.4g}\n".format(m,b),c='deepskyblue')
    print("Parameters STED resolution: \nalpha={} \ndc={}\nr2={}".format(m,b,rsquared))
except RuntimeError:
    print("Fit didn't work")
    
def monoExp(x, m, t, b):
    """
    Calculates the value of a monotonically decreasing exponential function at a given point.
    
    Parameters:
    x (float): The input value.
    m (float): The amplitude of the exponential function.
    t (float): The decay rate of the exponential function.
    b (float): The baseline value of the exponential function.
    
    Returns:
    float: The value of the monotonically decreasing exponential function at the given point.
    """
    return m * np.exp(-t * x) + b

xs = np.linspace(np.min(pp), np.max(pp),(np.max(pp)-np.min(pp)).astype(int))
# Fit the SPLIT-STED resolution vs STED power graph
try:
    params, cv = curve_fit(monoExp,data_objectifs['STEDpercent'],data_objectifs['Res splitsted']*20,p0=(dfmeanoverall['Res splitsted']*20,3.4E-3,dfmeanoverall['Res splitsted']*20),maxfev = 10000)
    print("params",params)
    m, t, b = params
    y_pred = monoExp(data_objectifs['STEDpercent'], m,t,b)
    rsquared=sklearn.metrics.r2_score(data_objectifs['Res splitsted']*20, y_pred)
    Ax.plot(xs, monoExp(xs,m,t,b),'-',c='mediumblue')
    Ax.text(0,125,"Parameters:\nA={0:.3g}\nt={1:.4g}\nb={2:.3g}".format(m,t,b),c='mediumblue')
    print("Parameters SPLIT Resolution:\nA={}\nt={}\nb={}\n r2={}".format(m,t,b,rsquared))
except RuntimeError:
    print("Fit didn't work")

# Plot the photobleaching vs STED power graph
Fig2,Ax2=plt.subplots(figsize=(8,6))
Ax2.errorbar(x=dfmean['STEDpercent'],y=dfmean['Bleach']*100,yerr=dfstd['Bleach']*100, fmt="o",c='gold',label='Photobleaching')
Ax2.set_ylabel('Photobleaching [%]',color='gold')
Ax.set_xlabel('STED depletion power [mW]')
Ax2.set_xlabel('STED depletion power [%]')
Ax.set_ylabel('Resolution [nm]')
Ax.legend(loc='upper right')
Ax2.legend(loc='upper left')
Ax.set_ylim([75,325])
mwpowers=["0","44","88","132","176"]
Ax.set_xticks([0,10,20,30,40])
Ax.set_xticklabels(mwpowers)
Ax.tick_params(axis='y', labelsize=16)

# Save the resolution vs STED power and photobleaching vs STED power graphs
Fig.savefig(os.path.join(savefolder,'PowervsResolution_SPLITSTED_Cy3_DTCWT_GroundTruth.png'),transparent='True', bbox_inches="tight")
Fig.savefig(os.path.join(savefolder,'PowervsResolution_SPLITSTED_Cy3_DTCWT_GroundTruth.pdf'),transparent='True', bbox_inches="tight")
Fig2.savefig(os.path.join(savefolder,'PowervsBleach_SPLITSTED_Cy3_DTCWT.png'),transparent='True', bbox_inches="tight")
Fig2.savefig(os.path.join(savefolder,'PowervsBleach_SPLITSTED_Cy3_DTCWT.pdf'),transparent='True', bbox_inches="tight")

plt.show()

