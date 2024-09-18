

""""   
    This script is used to compute the lifetime of the entire foreground in each FLIM  image in a folder
    The lifetime is measured using a fit of the histogram with a biexponential model with MLE error estimation
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
from Multi_fit import ExpFun_bi_MLE_tau_and_alpha
from Main_functions import get_foreground,load_msr
from lifetime import LifetimeOverlayer
import seaborn


# -----------------------------------------------------------

# Path to the folder containing the images

#filename =easygui.diropenbox(default=os.path.expanduser("~Desktop"))
filename= os.path.join('T:', os.sep,'adeschenes','SimulationDataset_STEDFLIM','Cy3',"PSD95_STORANGE")
#filename= os.path.join('T:', os.sep,'adeschenes','SimulationDataset_STEDFLIM','Cy3',"rabBassoon_CF594")

# Dictionary of the image identifiers (Channel names) to be included

#mapcomp = {'Conf635': 'Conf_635P {2}','STED635': 'STED_635P {2}'}
mapcomp = {'Confocal':'Confocal_561 {11}', 'STED':'STED 561 {11}'}

# Make list of all the images in the folder
print(filename)
path=os.path.join(filename,"*.msr")
images = glob.glob(path)
print('There are ',len(images), 'Images in this folder')

# Ask the user for a name for the output folder and create it
savefoldername=str(input("Name of Output folder: "))
savefolder = os.path.join(os.path.expanduser("~/Desktop"), "MLEForeground_"+savefoldername)
os.makedirs(savefolder, exist_ok=True)

# Create a dataframe to store the lifetime values
Overall_data = pd.DataFrame(
    columns=["Power", 'image',"tau1", "tau2", "f1", "f2","tau_mean"])

# Initial guess for the biexponential model
x0 = [3.16, 0.5, 0.5, 0.5] #tau1, tau2, f1, f2 Homer Orange
#x0 = [2.5, 0.5, 0.5, 0.5] #tau1, tau2, f1, f2 Bassoon CF594

# -----------------------------------------------------------
#     Choose the color for the graphs
graphcolor="hotpink"
#graphcolor='deepskyblue'
# -----------------------------------------------------------




for image_id,imagei in enumerate(images):
    print("##################")
    print("Image {} of {}".format(image_id,len(images)))
    print("##################")
        # Use the last part of the image name to get the STED power
    sted_percent = str(os.path.basename(imagei).split('_')[-1].split('percentSTED')[0])
    conf_percent=0
    print(os.path.basename(imagei))
    imagemsr = load_msr(imagei)

# -----------------------------------------------------------
#     Open mapcomp's images

    for k,key in enumerate(mapcomp):
        print(mapcomp[key])
        image1=imagemsr[mapcomp[key]]
        dim = image1.shape
        if k==0:
            ov_data = [int(conf_percent), os.path.basename(imagei)]
        else:
            ov_data=[int(sted_percent),os.path.basename(imagei)]

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
    #   Computes lifetime with MLE fit of biexponential model

        w = minimize(ExpFun_bi_MLE_tau_and_alpha, 
                    x0 = x0,
                    args = [absci, y], 
                    bounds = [(0,5), (0, 1),(0,1), (0,1)]
                    ) # (min, max) : tau1, tau2, f1, f2

        print(w)
        p1 = w['x'][0]  #tau1
        p2 = w['x'][1]  #tau2
        p3 = w['x'][2]  #f1
        p4 = w['x'][3]  #f2 
        print('tau1 :', p1)
        print('tau2 :', p2)
        print('f1 :', p3)
        print('f2 :', p4)
        # -----------------------------------

        # Add the lifetime values to the dataframe
        ov_data.extend([p1,p2,p3,p4, (p1 * p3) + (p2 * p4)])
        Overall_data.loc[len(Overall_data)] = ov_data
Overall_data.to_csv(os.path.join(savefolder, "MLE_foreground_Biexp_{}.csv".format(savefoldername)))

print(Overall_data.shape)
print(list(Overall_data.columns))

# -----------------------------------------------------------
#    Plot the lifetime values as a function of STED power
dfmeanconf=Overall_data.mean(numeric_only=True)
dfstdconf=Overall_data.std()

dfmean=Overall_data.groupby("Power", as_index=False).mean(numeric_only=True)
dfstd=Overall_data.groupby("Power", as_index=False).std()

print(dfmean.shape)
print(list(dfmean.columns))
print(Overall_data.shape)
print(list(Overall_data.columns))

Fig,Ax=plt.subplots(figsize=(6,4))
Fig1,Ax1=plt.subplots(figsize=(6,4))
Fig2,Ax2=plt.subplots(figsize=(6,4))
Fig3,Ax3=plt.subplots(figsize=(3,2))
Ax3.errorbar(x=dfmean["Power"],y=dfmean["tau_mean"],yerr=dfstd["tau_mean"], fmt="o",c=graphcolor,label='STED',ecolor='k',capsize=5,elinewidth=2)
#Ax.errorbar(x=0,y=dfmeanconf["lifetime conf"],yerr=dfstdconf["lifetime conf"], fmt="o",c= 'mediumblue',label='Confocal',ecolor='mediumblue',capsize=10,elinewidth=2,ms=15)
seaborn.violinplot(data=Overall_data,x=Overall_data["Power"],y=Overall_data["tau1"],ax=Ax,width=0.7)
seaborn.violinplot(data=Overall_data,x=Overall_data["Power"],y=Overall_data["tau2"],ax=Ax,width=0.7)
seaborn.violinplot(data=Overall_data,x=Overall_data["Power"],y=Overall_data["f1"],ax=Ax1,width=0.7)
seaborn.violinplot(data=Overall_data,x=Overall_data["Power"],y=Overall_data["f2"],ax=Ax1,width=0.7)
seaborn.violinplot(data=Overall_data,x=Overall_data["Power"],y=Overall_data["tau_mean"],ax=Ax2,width=0.7)
#seaborn.violinplot(data=Overall_data,x=0*Overall_data["Power"],y=Overall_data["lifetime conf"],ax=Ax,width=0.7)
Ax.set_ylim([0,4])
Ax1.set_ylim([0,1])
Ax2.set_ylim([0,4])
Ax3.set_ylim([0,4])

Fig.savefig(os.path.join(savefolder,'LifetimevsSTEDPower_{}.pdf'.format(savefoldername)), transparent='True', bbox_inches="tight")
Fig1.savefig(os.path.join(savefolder,'LifetimevsSTEDPower_Fractions_{}.pdf'.format(savefoldername)), transparent='True', bbox_inches="tight")
Fig2.savefig(os.path.join(savefolder,'ViolinPlot_LifetimevsSTEDPower_Mean_{}.pdf'.format(savefoldername)), transparent='True', bbox_inches="tight")


Fig3.savefig(os.path.join(savefolder,'LifetimevsSTEDPower_Mean_{}.pdf'.format(savefoldername)), transparent='True', bbox_inches="tight")

plt.show()
