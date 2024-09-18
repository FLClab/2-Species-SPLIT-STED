
"""
Program to read csv files containing SPLIT-STED metrics 
produced by: SPLIT-STED_CWF_multiSTEDpercent_Metrics.py

and 
Fluorescence lifetime measurements produced by histogram fitting with a bi-exponential model produced by:
Hist_biexp_fit_MLE_foreground_AllFolder.py

Create plots correlating lifetime values and differences and SPLIT-STED resolution metrics as a function of STED power
"""



import os.path
from sys import path as path1;

dossier = os.path.expanduser("~/Documents/Github/Abberior-STED-FLIM/Functions")
path1.append(dossier)
import functools
from statistics_functions import get_significance
import matplotlib.pyplot as plt
import numpy
import glob
import itertools
import tifffile
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
import seaborn
import pandas as pd
import scipy
from sklearn.linear_model import LinearRegression
import skimage
import matplotlib
matplotlib.rcParams['axes.linewidth'] = 0.8
# Labels for the different samples
labelsSPLIT=["Bassoon","Homer"]

# Folders containing the SPLIT-STED metrics for the 2 different samples
folder2=os.path.join(os.path.expanduser("Desktop"),"SPLIT_Homer_MediumAcq")
folder1=os.path.join(os.path.expanduser("Desktop"),"SPLIT_Bassoon_MediumAcq")

folder2=os.path.join(os.path.expanduser("Desktop"),"PSD95_Orange_MediumPlus_SPLIT")
folder1=os.path.join(os.path.expanduser("Desktop"),"rabBassoon_CF594_MediumPlus_SPLIT")
folder1=os.path.join("D:",os.sep,"FLIM_MediumAcq","PSD-Bassoon_Cy3","SPLIT-STED","rabBassoon_CF594_MediumPlus_SPLIT")
folder2=os.path.join("D:",os.sep,"FLIM_MediumAcq","PSD-Bassoon_Cy3","SPLIT-STED","PSD95_Orange_MediumPlus_SPLIT")
foldersSPLIT=[folder1,folder2]

# Folders containing the lifetime measurements for the 2 different samples
folder1=os.path.join(os.path.expanduser("Desktop"),"MLEForeground_Bassoon_CF594_BiExp")
folder2=os.path.join(os.path.expanduser("Desktop"),"MLEForeground_Homer_Orange_BiExp")
folder1=os.path.join(os.path.expanduser("Desktop"),"MLEForeground_rabBassoon_CF594_BiExp")
folder2=os.path.join(os.path.expanduser("Desktop"),"MLEForeground_PSD95_STORANGE_BiExp")

foldersMLE=[folder1,folder2]

# Create the figures and axes for the plots
Fig, Ax = plt.subplots(nrows=2,ncols=2,figsize=(6, 6),sharex=True)
Fig3, Ax3 = plt.subplots(figsize=(6,4),sharex=True)
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
colors=[["#f9a3cbff","#ef87beff","#e569b3ff","#bf4290ff",'magenta'],["#55d0ffff","#00acdfff","#0080bfff","#00456bff",'cyan']]

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


# Read the csv files containing the lifetime measurements and combine them into a single dataframe
cumdf=[]
csvfull=[]
for i,folder in enumerate(foldersMLE):

    #colors = ["xkcd:peacock blue", "xkcd:brick orange"]
    csvlist=glob.glob(os.path.join(folder,"*.csv"))


    Overall_data= pd.concat(map(pd.read_csv, csvlist))
    cumdf.append(Overall_data)

    print(Overall_data.shape)
    print(list(Overall_data.columns))
    
Overall_data_MLE= pd.concat([df.assign(identity=k) for k,df in zip(labelsSPLIT,cumdf)])
print(Overall_data_MLE.shape)
MLEcumulative_Mean=[]

# Plot the resolution of STED and SPLIT-STED images as a function of STED power for each sample
Ax[0,0].scatter(Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[1]]["STEDpercent"],Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[1]]['Res sted stack']*20)    
Ax[0,0].scatter(Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[1]]["STEDpercent"],Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[1]]['Res splitsted']*20)  
Ax[0,1].scatter(Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[0]]["STEDpercent"],Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[0]]['Res sted stack']*20)    
Ax[0,1].scatter(Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[0]]["STEDpercent"],Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[0]]['Res splitsted']*20) 

# Plot the lifetime values as a function of STED power for each sample
Ax[1,0].scatter(Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[1]]["Power"],Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[1]][ 'tau_mean'])    
#Ax[1,1].scatter(Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[1]]["Power"],Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[1]]['f2'])  
Ax[1,0].scatter(Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[0]]["Power"],Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[0]][ 'tau_mean'])    

#Ax[1,1].scatter(Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[0]]["Power"],Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[0]]['f2']) 

powers = numpy.unique(Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[1]]["Power"])
table_deltatau=[]
table_deltatau_median=[]
table_deltatau_std=[]
table_deltatau_temp=[]

# Calculate the difference in lifetime values between all possible combinations of images of the 2 samples for each STED power. 
#Calculate the median and standard deviation of these differences
for power in powers:
    table_deltatau_temp=[]
    pairs = list(itertools.product(Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0]) & (Overall_data_MLE["Power"] == power)][ 'tau_mean'], Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1]) &(Overall_data_MLE["Power"] == power)][ 'tau_mean']))
    print(Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0]) & (Overall_data_MLE["Power"] == power)][ 'tau_mean'].shape)
    print(Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1]) & (Overall_data_MLE["Power"] == power)][ 'tau_mean'].shape)
    print(len(pairs))
    for pair in pairs:
        table_deltatau_temp.append([power,numpy.abs(pair[0]-pair[1])])
    table_deltatau.extend(table_deltatau_temp)
    table_deltatau_median.append(numpy.median(table_deltatau_temp,axis=0))
    table_deltatau_std.append(numpy.std(table_deltatau_temp,axis=0))

table_deltatau=numpy.array(table_deltatau)
table_deltatau_median=numpy.array(table_deltatau_median)
table_deltatau_std=numpy.array(table_deltatau_std)
print(table_deltatau_median)
print(table_deltatau_std)

# Plot the difference in lifetime values between the 2 samples as a function of STED power

#Ax[1,1].scatter(table_deltatau[:,0],table_deltatau[:,1]) 
seaborn.violinplot(x=table_deltatau[:,0],y=table_deltatau[:,1],ax=Ax3,width=0.7)
Fig2, Ax2 = plt.subplots(nrows=2,ncols=3,figsize=(12,8))
powers=[10,20,30,40]
#powers = [5,10,15,20]


tauerror=table_deltatau_std[table_deltatau_median[:,0]==0][0,1]
tau=table_deltatau_median[table_deltatau_median[:,0]==0][0,1]

# Calculate the mean and standard deviation of the resolution metrics for the SPLIT-STED images as a function of STED power for each sample
Mean_Power_SPLIT_list=[]
STD_Power_SPLIT_list=[]
for i,id in enumerate(labelsSPLIT):
    Mean_Power_SPLIT=Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == id].groupby(by=["STEDpercent"]).mean(numeric_only=True).reset_index()
    STD_Power_SPLIT=Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == id].groupby(by=["STEDpercent"]).std().reset_index()
    Mean_Power_SPLIT_list.append(Mean_Power_SPLIT)
    STD_Power_SPLIT_list.append(STD_Power_SPLIT)

# Plot vertical lines indicating the lifetime difference between the 2 samples for confocal images
    Ax2[i,1].axvline(tau,label="Confocal",c="orange")
    Ax2[i,1].axvspan(tau-tauerror,tau+tauerror,alpha=0.3,facecolor="orange")
    

    for p,power in enumerate(powers):
        print(table_deltatau_std[table_deltatau_median[:,0]==power][0,1])
        xerror=table_deltatau_std[table_deltatau_median[:,0]==power][0,1]
        
        a=Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res sted stack'].item()*20

        b=STD_Power_SPLIT.loc[STD_Power_SPLIT["STEDpercent"]==power]['Res sted stack'].item()*20
        print(a,b)

        Ax2[i,1].axhspan(a-b,a+b,alpha=0.3,facecolor=colors[i][p])
        Ax2[i,1].axhline(a,label="pSTED={}%".format(power),c=colors[i][p])


        #Ax2[i,0].errorbar(x=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],xerr=xerror,y=Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res sted stack']*20,yerr=STD_Power_SPLIT.loc[STD_Power_SPLIT["STEDpercent"]==power]['Res sted stack']*20,fmt="o",
        #           capsize=3.5, elinewidth=0.8, ms=10,ecolor="k",markeredgecolor="k",c=colors[i][p],label="pSTED={}%".format(power))
        Ax2[i,1].errorbar(x=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],xerr=xerror,y=Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res splitsted']*20,yerr=STD_Power_SPLIT.loc[STD_Power_SPLIT["STEDpercent"]==power]['Res splitsted']*20,fmt="o",
                   capsize=3.5, elinewidth=0.8, ms=10,ecolor="k",markeredgecolor="k",c=colors[i][p],label="pSTED={}%".format(power))
        Ax2[i,2].errorbar(x=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],xerr=xerror,y=(Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res sted stack']*20-Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res splitsted']*20)/Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res sted stack']*20,yerr=STD_Power_SPLIT.loc[STD_Power_SPLIT["STEDpercent"]==power]['Res splitsted']*20,fmt="o",
                   capsize=3.5, elinewidth=0.8, ms=10,ecolor="k",markeredgecolor="k",c=colors[i][p],label="pSTED={}%".format(power))


# Plot the confocal resolution of the 2 samples with the color of the circle indicating the lifetime difference between the 2 samples
Ax2[0,0].errorbar(x=Mean_Power_SPLIT_list[0]['Res conf stack'].mean()*20,xerr=STD_Power_SPLIT_list[0]['Res conf stack'].mean()*20,
                      y=Mean_Power_SPLIT_list[1]['Res conf stack'].mean()*20,yerr=STD_Power_SPLIT_list[1]['Res conf stack'].mean()*20,
                      fmt="o",capsize=3.5, elinewidth=0.8,ecolor="k",markeredgecolor="k",c="w",zorder=0)        
Ax2[0,0].scatter(x=Mean_Power_SPLIT_list[0]['Res conf stack'].mean()*20,
                      y=Mean_Power_SPLIT_list[1]['Res conf stack'].mean()*20,
                     c=table_deltatau_median[table_deltatau_median[:,0]==0][0,1],edgecolor="k",
                    s=100,vmin=0.75,vmax=tau,cmap="turbo",label="pSTED={}%".format(0),zorder=100) 
# Plot the STED resolution of the 2 samples with the color of the circle indicating the lifetime difference between the 2 samples
for p,power in enumerate(powers):
    Ax2[0,0].errorbar(x=Mean_Power_SPLIT_list[0].loc[Mean_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,xerr=STD_Power_SPLIT_list[0].loc[STD_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,
                      y=Mean_Power_SPLIT_list[1].loc[Mean_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20,yerr=STD_Power_SPLIT_list[1].loc[STD_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20,
                      fmt="o",capsize=3.5, elinewidth=0.8,ecolor="k",markeredgecolor="k",c="w",zorder=0)        
    scat=Ax2[0,0].scatter(x=Mean_Power_SPLIT_list[0].loc[Mean_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,
                     y=Mean_Power_SPLIT_list[1].loc[Mean_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20,
                     c=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],edgecolor="k",
                    s=100,vmin=0.75,vmax=tau,cmap="turbo",label="pSTED={}%".format(power),zorder=100) 
Fig2.colorbar(scat, ax=Ax2[0,0], location='top')
Ax2[0,0].legend()
Ax2[0,1].legend()
Ax2[1,0].legend()
Ax2[1,1].legend()
Ax3.set_ylim([0,4])
Fig2.savefig('Resolution_vs_MLEbiexp.pdf', transparent='True', bbox_inches="tight")
Fig3.savefig('MLEBiexp_DeltaTau.pdf', transparent='True', bbox_inches="tight")

plt.show()

