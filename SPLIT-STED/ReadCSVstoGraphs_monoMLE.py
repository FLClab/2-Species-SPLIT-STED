"""
Program to read csv files containing SPLIT-STED metrics 
produced by: SPLIT-STED_CWF_multiSTEDpercent_Metrics.py

and 
Fluorescence lifetime measurements produced by histogram fitting with a mono-exponential model produced by:
Hist_monoexp_MLE_Foreground_AllFolder.py

Create plots correlating lifetime values and differences and SPLIT-STED resolution metrics as a function of STED power
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
import matplotlib
matplotlib.rcParams['axes.linewidth'] = 0.8

# Labels for the different samples
labelsSPLIT=["Bassoon","PSD95"]



# Folders containing the lifetime measurements for the 2 different samples
folder1=os.path.join("C:",os.sep,"Users","ANDES399","Desktop","MLEForeground_BassoonCF594")
folder2=os.path.join("C:",os.sep,"Users","ANDES399","Desktop","MLEForeground_PSD95_Orange")


foldersMLE=[folder1,folder2]
# Create the plots
Fig, Ax = plt.subplots(nrows=2,ncols=2,figsize=(6, 6))
Fig3, Ax3 = plt.subplots(figsize=(6,4))
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
colors=[["#f9a3cbff","#ef87beff","#e569b3ff","#bf4290ff"],["#55d0ffff","#00acdfff","#0080bfff","#00456bff"]]

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

# Plot the lifetime values as a function of STED power for each sample

Ax[1,0].scatter(Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[1]]["Power"],Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[1]][ 'lifetime'])    
#Ax[1,1].scatter(Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[1]]["Power"],Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[1]]['f2'])  
Ax[1,0].scatter(Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[0]]["Power"],Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[0]][ 'lifetime'])    

#Ax[1,1].scatter(Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[0]]["Power"],Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[0]]['f2']) 

powers = numpy.unique(Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[1]]["Power"])
table_deltatau=[]
table_deltatau_mean=[]
table_deltatau_std=[]
table_deltatau_temp=[]
# Calculate the difference in lifetime values between all possible combinations of images of the 2 samples for each STED power. 
#Calculate the median and standard deviation of these differences
for power in powers:
    table_deltatau_temp=[]
    pairs = list(itertools.product(Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0]) & (Overall_data_MLE["Power"] == power)][ 'lifetime'], Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1]) &(Overall_data_MLE["Power"] == power)][ 'lifetime']))
    print(Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0]) & (Overall_data_MLE["Power"] == power)][ 'lifetime'].shape)
    print(Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1]) & (Overall_data_MLE["Power"] == power)][ 'lifetime'].shape)
    print(len(pairs))
    for pair in pairs:
        table_deltatau_temp.append([power,numpy.abs(pair[0]-pair[1])])
    table_deltatau.extend(table_deltatau_temp)
    table_deltatau_mean.append(numpy.mean(table_deltatau_temp,axis=0))
    table_deltatau_std.append(numpy.std(table_deltatau_temp,axis=0))

table_deltatau=numpy.array(table_deltatau)
table_deltatau_mean=numpy.array(table_deltatau_mean)
table_deltatau_std=numpy.array(table_deltatau_std)
print(table_deltatau_mean)
print(table_deltatau_std)
# Plot the difference in lifetime values between the 2 samples as a function of STED power
seaborn.violinplot(x=table_deltatau[:,0],y=table_deltatau[:,1],ax=Ax3,width=0.7)
tauerror=table_deltatau_std[table_deltatau_mean[:,0]==0][0,1]
tau=table_deltatau_mean[table_deltatau_mean[:,0]==0][0,1]
Fig4, Ax4 = plt.subplots(nrows=2,ncols=3,figsize=(12,8))
Fig2, Ax2 = plt.subplots(figsize=(4,3))
Fig5, Ax5 = plt.subplots(figsize=(3,2))
powers=[10,20,30,40]

Fig3.savefig('MLEmonoexp_DeltaTau.pdf', transparent='True', bbox_inches="tight")
Fig5.savefig('MLEmonoexp_TauBoxplot.pdf', transparent='True', bbox_inches="tight")

# Plot the lifetime values as a function of STED power for each sample
seaborn.boxplot(x=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1])]["Power"], y=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1])]["lifetime"],ax=Ax5, showfliers = False, whis=(0, 100),boxprops={"facecolor": None, "edgecolor":"hotpink"},medianprops={"color": "k"})
seaborn.boxplot(x=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0])]["Power"], y=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0])]["lifetime"],ax=Ax5, showfliers = False, whis=(0, 100),boxprops={"facecolor": None, "edgecolor":"deepskyblue"},medianprops={"color": "k"})


plt.show()