"""
Program to read csv files containing SPLIT-STED metrics 
produced by: SPLIT-STED_CWF_multiSTEDpercent_Metrics.py

and 
Phasor Ellipses produced by:
PhasorDistribution_STED_EllipseFit_2Species.py

Create plots correlating distance between ellipses and SPLIT-STED resolution metrics as a function of STED power
"""


import os.path
from sys import path as path1;
dossier = os.path.expanduser("~/Documents/Github/2-Species-SPLIT-STED/Functions")
path1.append(dossier)
import functools
from statistics_functions import get_significance
import matplotlib.pyplot as plt
import numpy
import glob
import itertools
import tifffile
import seaborn
import pandas as pd
import scipy
import skimage
import matplotlib
matplotlib.rcParams['axes.linewidth'] = 0.8
# Labels for the different samples
labelsSPLIT=["Bassoon","Homer"]
# Folders containing the SPLIT-STED metrics for the 2 different samples (1 csv file per sample)
folder1=os.path.join(os.path.expanduser("Desktop"),"SPLIT_Homer_MediumAcq")
folder2=os.path.join(os.path.expanduser("Desktop"),"SPLIT_Bassoon_MediumAcq")
folder2=os.path.join(os.path.expanduser("Desktop"),"PSD95_Orange_MediumPlus_SPLIT")
folder1=os.path.join(os.path.expanduser("Desktop"),"rabBassoon_CF594_MediumPlus_SPLIT")
#folder1=os.path.join(os.path.expanduser("Desktop"),"Bassoon_ST635P_SPLIT")
#folder2=os.path.join(os.path.expanduser("Desktop"),"PSD95_AF647_SPLIT")
foldersSPLIT=[folder1,folder2]

# Folder containing the phasor ellipse data (1 csv file for pair of samples)
folderellipse=os.path.join(os.path.expanduser("Desktop"),"BassoonCF594_HomerOrange_MediumAcq_PhasorDists")
folderellipse=os.path.join(os.path.expanduser("Desktop"),"rabBassoonCF594_PSD95Orange_PhasorDists")

# Create the figures and axes for the plots
Fig, Ax = plt.subplots(nrows=2,ncols=2,figsize=(6, 6),sharex=True,sharey='row')
Ax[0,0].set_ylim([45, 320])
Ax[0,1].set_ylim([45, 320])
mwpowers=["0","44","88","132","176"]
Ax[0,0].tick_params(axis ='both',length=2, width=0.8)
Ax[0,1].tick_params(axis ='both',length=2, width=0.8)
Ax[1,0].tick_params(axis ='both',length=2, width=0.8)
Ax[1,1].tick_params(axis ='both',length=2, width=0.8)
Ax[0,1].set_xticks([0,10,20,30,40])
Ax[0,1].set_xticklabels(mwpowers,fontsize=16)
cumdf=[]
csvfull=[]
colors = ["xkcd:peacock blue", "xkcd:brick orange"]
colors=[["#f9a3cbff","#ef87beff","#e569b3ff","#bf4290ff"],["#55d0ffff","#00acdfff","#0080bfff","#00456bff"]]

# Read the csv files containing the SPLIT-STED metrics and combine them into a single dataframe
for i,folder in enumerate(foldersSPLIT):

    
    csvlist=glob.glob(os.path.join(folder,"*.csv"))
    Overall_data= pd.concat(map(pd.read_csv, csvlist))
    cumdf.append(Overall_data)
    print(Overall_data.shape)
    print(list(Overall_data.columns))
    
Overall_data_SPLIT= pd.concat([df.assign(identity=k) for k,df in zip(labelsSPLIT,cumdf)])
print(Overall_data_SPLIT.shape)
print(Overall_data_SPLIT.columns)

cumdf=[]
csvfull=[]
# Read the csv files containing the ellipse data and combine them into a single dataframe
csvlist=glob.glob(os.path.join(folderellipse,"*.csv"))
Overall_data_ellipse= pd.concat(map(pd.read_csv, csvlist))
print(Overall_data_ellipse)
   
    
MLEcumulative_Mean=[]

# Plot the resolution of STED and SPLIT-STED images as a function of STED power for each sample
Ax[0,0].scatter(Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[1]]["STEDpercent"],Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[1]]['Res sted stack']*20)    
Ax[0,0].scatter(Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[1]]["STEDpercent"],Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[1]]['Res splitsted']*20)  
Ax[0,1].scatter(Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[0]]["STEDpercent"],Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[0]]['Res sted stack']*20)    
Ax[0,1].scatter(Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[0]]["STEDpercent"],Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[0]]['Res splitsted']*20) 

#Plot the shortest distance between the ellipses as a function of STED power
Ax[1,0].scatter(Overall_data_ellipse["Power"],Overall_data_ellipse[ 'Shortest Distance'])    
#Plot the intersection over union between the ellipses as a function of STED power
Ax[1,1].scatter(Overall_data_ellipse["Power"],Overall_data_ellipse[ 'IOU'])  


powers = numpy.unique(Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[1]]["STEDpercent"])


# Create the figures and axes for the plots
Fig2, Ax2 = plt.subplots(nrows=2,ncols=3,figsize=(6, 6),sharex=True,sharey='col')
powers=[10,20,30,40]
#powers=[5,10,15,20]

# Calculate the mean and standard deviation of the ellipse metrics for each STED power
Mean_Power_SPLIT_list=[]
STD_Power_SPLIT_list=[]
Mean_Power_Ellipse=Overall_data_ellipse.groupby(by=["Power"]).mean(numeric_only=True).reset_index()
STD_Power_Ellipse=Overall_data_ellipse.groupby(by=["Power"]).std().reset_index()


# Calculate the mean and standard deviation of the SPLIT-STED metrics for each sample and STED power
for i,id in enumerate(labelsSPLIT):
    Mean_Power_SPLIT=Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == id].groupby(by=["STEDpercent"]).mean(numeric_only=True).reset_index()
    STD_Power_SPLIT=Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == id].groupby(by=["STEDpercent"]).std().reset_index()
    Mean_Power_SPLIT_list.append(Mean_Power_SPLIT)
    STD_Power_SPLIT_list.append(STD_Power_SPLIT)
# Calculate the mean and standard deviation of the ellipse metrics for confocal and display as vertical rectangle on the plot
    ell=Mean_Power_Ellipse.loc[Mean_Power_Ellipse["Power"]==0]['Shortest Distance'].item()
    ellerr=STD_Power_Ellipse.loc[Mean_Power_Ellipse["Power"]==0]['Shortest Distance'].item()
    Ax2[i,1].axvline(ell,label="Confocal",c="orange")
    Ax2[i,1].axvspan(ell-ellerr,ell+ellerr,alpha=0.3,facecolor="orange")
# Graph the STED resolution of a sample vs the ellipse distance between the two samples
    for p,power in enumerate(powers):
        print(Mean_Power_Ellipse.loc[Mean_Power_Ellipse["Power"]==power]['Shortest Distance'])


        a=Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res sted stack'].item()*20

        b=STD_Power_SPLIT.loc[STD_Power_SPLIT["STEDpercent"]==power]['Res sted stack'].item()*20

        
        Ax2[i,1].axhspan(a-b,a+b,alpha=0.3,facecolor=colors[i][p])
        Ax2[i,1].axhline(a,label="pSTED={}%".format(power),c=colors[i][p])

        Ax2[i,0].errorbar(x=Mean_Power_Ellipse.loc[Mean_Power_Ellipse["Power"]==power]['Shortest Distance'],xerr=STD_Power_Ellipse.loc[Mean_Power_Ellipse["Power"]==power]['Shortest Distance'],y=Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res sted stack']*20,yerr=STD_Power_SPLIT.loc[STD_Power_SPLIT["STEDpercent"]==power]['Res sted stack']*20,fmt="o",
                   capsize=3.5, elinewidth=0.8, ms=10,ecolor="k",markeredgecolor="k",c=colors[i][p],label="pSTED={}%".format(power))
        Ax2[i,1].errorbar(x=Mean_Power_Ellipse.loc[Mean_Power_Ellipse["Power"]==power]['Shortest Distance'],xerr=STD_Power_Ellipse.loc[Mean_Power_Ellipse["Power"]==power]['Shortest Distance'],y=Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res splitsted']*20,yerr=STD_Power_SPLIT.loc[STD_Power_SPLIT["STEDpercent"]==power]['Res splitsted']*20,fmt="o",
                   capsize=3.5, elinewidth=0.8, ms=10,ecolor="k",markeredgecolor="k",c=colors[i][p],label="pSTED={}%".format(power))
        Ax2[i,2].errorbar(x=Mean_Power_Ellipse.loc[Mean_Power_Ellipse["Power"]==power]['Shortest Distance'],xerr=STD_Power_Ellipse.loc[Mean_Power_Ellipse["Power"]==power]['Shortest Distance'],y=(Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res sted stack']*20-Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res splitsted']*20)/Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res sted stack']*20,yerr=STD_Power_SPLIT.loc[STD_Power_SPLIT["STEDpercent"]==power]['Res splitsted']*20,fmt="o",
                   capsize=3.5, elinewidth=0.8, ms=10,ecolor="k",markeredgecolor="k",c=colors[i][p],label="pSTED={}%".format(power))
# Graph the STED resolution of both samples with the ellipse distance as color      
for p,power in enumerate(powers):
    Ax2[0,0].errorbar(x=Mean_Power_SPLIT_list[0].loc[Mean_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,xerr=STD_Power_SPLIT_list[0].loc[STD_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,
                      y=Mean_Power_SPLIT_list[1].loc[Mean_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20,yerr=STD_Power_SPLIT_list[1].loc[STD_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20,
                      fmt="o",capsize=3.5, elinewidth=0.8,ecolor="k",markeredgecolor="k",c="w",zorder=0)        
    Ax2[0,0].scatter(x=Mean_Power_SPLIT_list[0].loc[Mean_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,
                     y=Mean_Power_SPLIT_list[1].loc[Mean_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20,
                     c=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],edgecolor="k",
                    s=100,vmin=0.6,vmax=1,cmap="jet",label="pSTED={}%".format(power),zorder=100) 
Ax2[0,0].legend()
Ax2[0,1].legend()
Ax2[1,0].legend()
Ax2[1,1].legend()
Ax2[0,2].legend()
Ax2[1,2].legend()

plt.show()

