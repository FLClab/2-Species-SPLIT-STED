
""""   
    This script is used to compute the lifetime of the entire foreground in each FLIM  image in a folder
    The lifetime is measured using a fit of the histogram with a mono-exponential model with MLE error estimation
The script will output a csv file with the lifetime values for each image in the folder
"""

import os
import matplotlib.pyplot as plt
import glob
import numpy
import easygui

from scipy.optimize import minimize
from tqdm import tqdm
import tifffile
import os.path
import pandas as pd
from sys import path as path1;
dossier = os.path.expanduser("~/Documents/Github/2-Species-SPLIT-STED/Functions")
path1.append(dossier)
from Mono_fit import ExpFunMono_MLE
from Main_functions import get_foreground,load_image,select_channel

from tiffwrapper import LifetimeOverlayer
import seaborn
# -----------------------------------------------------------


graphcolor="deepskyblue"
#graphcolor="hotpink"


# Path to the folder containing the images
#filename =easygui.diropenbox(default=os.path.expanduser("~Desktop"))

filename= os.path.join('T:', os.sep,'adeschenes','SimulationDataset_STEDFLIM','Cy3',"PSD95_STORANGE")
filename= os.path.join('T:', os.sep,'adeschenes','SimulationDataset_STEDFLIM','Cy3',"rabBassoon_CF594")
#filename= os.path.join('T:', os.sep,'adeschenes','SimulationDataset_STEDFLIM','Cy5','alphaTubulin_Alexa647')
#filename= os.path.join('T:', os.sep,'adeschenes','SimulationDataset_STEDFLIM','Cy5','B2Spectrin Alexa647')
#filename= os.path.join('T:', os.sep,'adeschenes','SimulationDataset_STEDFLIM','Cy5','rabBassoon STAR635P')
#filename=  os.path.join('U:', os.sep,'adeschenes',"2023-12-21_FLIM_MediumAcq_Spectrin_Actin_Bassoon","Bassoon_CF594_STEDPowerBleach_MediumAcq_1")
#filename= os.path.join('U:', os.sep,'adeschenes',"2023-12-21_FLIM_MediumAcq_Spectrin_Actin_Bassoon","B2Spectrin_STOrange_STEDPowerBleach_MediumAcq_1")
#filename=os.path.join('T:', os.sep,'adeschenes',"SimulationDataset_STEDFLIM","Cy3","Bassoon_CF594","HighP")
#filename=os.path.join('T:', os.sep,'adeschenes',"SimulationDataset_STEDFLIM","Cy3","Homer_STORANGE","HighP")
# Dictionary of the image identifiers (Channel names) to be included
mapcomp = {'CONF561': 'Confocal_561 {11}', 'STED561' : 'STED 561 {11}'}
mapcomp = {'CONF561': 0, 'STED561' : 1}
#mapcomp = {'Conf635': 'Conf_635P {2}','STED635': 'STED_635P {2}'}
filename =easygui.diropenbox(default=os.path.expanduser("~Desktop"))

# Make list of all the images in the folder
extension = ".msr"
path=os.path.join(filename,"*.msr")
images = glob.glob(path)
print('There are ',len(images), ' msr files in this folder')
if len(images) == 0:
    path=os.path.join(filename,"*.tiff")
    images = glob.glob(path)
    print('There are ',len(images), ' tiff files in this folder')
    extension = ".tiff"

# Ask the user for a name for the output folder and create it
savefoldername=str(input("Name of Output folder: "))
savefolder = os.path.join(os.path.expanduser("~/Desktop"), "MLEForeground_"+savefoldername)
os.makedirs(savefolder, exist_ok=True)

# Create a dataframe to store the lifetime values
Overall_data = pd.DataFrame(
    columns=["Power", 'image',"lifetime"])

# -----------------------------------------------------------




for image_id,imagei in enumerate(images):
    print("##################")
    print("Image {} of {}\n".format(image_id,len(images)))
    print("##################")
    # Use the last part of the image name to get the STED power
    sted_percent = str(os.path.basename(imagei).split('_')[-1].split('percentSTED')[0])
    conf_percent=0
    ov_data=[int(sted_percent),os.path.basename(imagei)]
    ov_data_conf = [int(conf_percent), os.path.basename(imagei)]
    print(os.path.basename(imagei))
    imagemsr = load_image(imagei)

# -----------------------------------------------------------
#     Open mapcomp's images

    for k,key in enumerate(mapcomp):
        print(mapcomp[key])
        image1=select_channel(imagemsr, mapcomp[key])
        #image1=imagemsr[mapcomp[key]]
        dim = image1.shape


        if dim[2] > 250 :
            image1 = image1[:,:,:250].astype(numpy.int16)
            print(image1.shape)
            dim=image1.shape

        dim = image1.shape

        imsum= numpy.sum(image1, axis=2, dtype = numpy.int16)

    # -----------------------------------------------------------

        #seuil = get_foreground(image1)
        seuil= 3
    # Sum of all the histograms of the foreground pixels
        y=numpy.sum(image1[imsum>seuil, :],axis=0)
    # Cut histogram to start at max value
        maxy = numpy.max(y)
        indice = numpy.argmax(y)
        y = y[indice:]
        y= y / y.sum()
        absci = numpy.linspace(0,y.shape[0]-1, num =y.shape[0])*0.08
#   Computes lifetime for sum of foreground with MLE fit of monoexponential model

        tau = 2
        nb_photons = 100
        bounds = [(0,0.0001),(0.1,5),(0,900000)] # (min, max): Background, tau, amp
        w = minimize(ExpFunMono_MLE,
                    x0 = [0,tau,nb_photons],
                    args = [absci, y],
                    bounds = bounds)
        bg, tau, amp =  w["x"]
        lifetime= w['x'][1]
        print(lifetime)
        if k==0:
            ov_data_conf.append(lifetime)
            Overall_data.loc[len(Overall_data)] = ov_data_conf
        else:
            ov_data.append(lifetime)
# Add the lifetime values to the dataframe
            Overall_data.loc[len(Overall_data)] = ov_data
Overall_data.to_csv(os.path.join(savefolder, "MLE_foreground_{}.csv".format(savefoldername)))

print(Overall_data.shape)
print(list(Overall_data.columns))

# -----------------------------------------------------------
#    Plot the lifetime values as a function of STED power
dfmeanconf=Overall_data.mean(numeric_only=True)
dfstdconf=Overall_data.std(numeric_only=True)

dfmean=Overall_data.groupby("Power", as_index=False).mean(numeric_only=True)
dfstd=Overall_data.groupby("Power", as_index=False).std(numeric_only=True)

print(dfmean.shape)
print(list(dfmean.columns))
print(Overall_data.shape)
print(list(Overall_data.columns))

Fig,Ax=plt.subplots(figsize=(3,2))
Fig2,Ax2=plt.subplots(figsize=(6,4))
#Ax.errorbar(x=dfmean["Power"],y=dfmean["lifetime"],yerr=dfstd["lifetime"], fmt="o",c=graphcolor,label='STED',ecolor='k',capsize=5,elinewidth=2)
#Ax.errorbar(x=0,y=dfmeanconf["lifetime conf"],yerr=dfstdconf["lifetime conf"], fmt="o",c= 'mediumblue',label='Confocal',ecolor='mediumblue',capsize=10,elinewidth=2,ms=15)
seaborn.violinplot(data=Overall_data,x=Overall_data["Power"],y=Overall_data["lifetime"],ax=Ax2,width=0.7)
#seaborn.violinplot(data=Overall_data,x=0*Overall_data["Power"],y=Overall_data["lifetime conf"],ax=Ax,width=0.7)

# Plot the lifetime for each depletion power with boxplots
seaborn.boxplot(
    data=Overall_data, x=Overall_data["Power"],y=Overall_data["lifetime"],ax=Ax)

# Add in points to show each observation
seaborn.stripplot( data=Overall_data, x=Overall_data["Power"],y=Overall_data["lifetime"],ax=Ax, size=4)

#Ax.set_ylim([1.5,3.5])
Ax.set_ylim([0,4])
Ax2.set_ylim([0,4])
Fig.savefig(os.path.join(savefolder,'LifetimevsSTEDPower_{}.pdf'.format(savefoldername)), transparent='True', bbox_inches="tight")
Fig2.savefig(os.path.join(savefolder,'ViolinPlot_LifetimevsSTEDPower_{}.pdf'.format(savefoldername)), transparent='True', bbox_inches="tight")
plt.show()