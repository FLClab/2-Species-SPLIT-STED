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
from statistics_functions import get_significance
import matplotlib.pyplot as plt
import numpy
import glob
import itertools
import tifffile
import seaborn
import pandas as pd
import scipy
import matplotlib
import easygui
matplotlib.rcParams['axes.linewidth'] = 0.8

# Labels for the different samples
labelsSPLIT=["Bassoon","PSD95"]
#labelsSPLIT=["B2Spectrin","Bassoon"]
#labelsSPLIT=["AlphaTubulin","Bassoon"]

# Folders containing the SPLIT-STED metrics for the 2 different samples
folder1=os.path.join(os.path.expanduser("Desktop"),"rabBassoon_CF594_SPLIT")
folder2=os.path.join(os.path.expanduser("Desktop"),"PSD95_STOrange_SPLIT")
foldersSPLIT=[folder1,folder2]

# Folders containing the lifetime measurements for the 2 different samples
folder1=os.path.join(os.path.expanduser("Desktop"),"MLEForeground_rabBassoon_CF594")
folder2=os.path.join(os.path.expanduser("Desktop"),"MLEForeground_PSD95_STORANGE")
foldersMLE=[folder1,folder2]

# Path of the folder where the results will be saved
savefolder=str(input("Name of Output folder: "))
savefolder = os.path.join(os.path.expanduser("~/Desktop"), savefolder)
os.makedirs(savefolder, exist_ok=True)

# Create the plots
Fig, Ax = plt.subplots(nrows=2,ncols=2,figsize=(6, 6))
Fig0, Ax0 = plt.subplots(nrows=1,ncols=2,figsize=(8, 4),sharey=True)
Fig3, Ax3 = plt.subplots(figsize=(6,4))
Ax[0,0].set_ylim([45, 320])
Ax[0,1].set_ylim([45, 320])
plt.setp(Ax0, ylim=(0,100))
#Ax0[0].set_ylim([0,100])
#Ax0[1].set_ylim([0,100])
mwpowers=["0","44","88","132","176"]
Ax[0,0].tick_params(axis ='both',length=2, width=0.8)
Ax[0,1].tick_params(axis ='both',length=2, width=0.8)
Ax[1,0].tick_params(axis ='both',length=2, width=0.8)
Ax[1,1].tick_params(axis ='both',length=2, width=0.8)
Ax[0,1].set_xticks([0,10,20,30,40])
Ax[0,1].set_xticklabels(mwpowers,fontsize=16)
Ax0[0].set_xticks([5,10,15,20,30])
Ax0[1].set_xticks([5,10,15,20,30])
mwpowers=["22","44","66","88","132"]
Ax0[0].set_xticklabels(mwpowers)
Ax0[1].set_xticklabels(mwpowers)
cumdf=[]
csvfull=[]
colors=[["#fcbcd7ff","#f9a3cbff","#ef87beff","#e569b3ff","#bf4290ff"],["#7ce8ffff","#55d0ffff","#00acdfff","#0080bfff","#00456bff"]]
#colors=[["#f9a3cbff","#ef87beff","#e569b3ff","#bf4290ff"],["#55d0ffff","#00acdfff","#0080bfff","#00456bff"]]
# Read the csv files containing the SPLIT-STED metrics and combine them into a single dataframe
for i,folder in enumerate(foldersSPLIT):

    #colors = ["xkcd:peacock blue", "xkcd:brick orange"]
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
#Ax0[0].scatter(Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[0]]["STEDpercent"],Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[0]]['Bleach']*100) 
#Ax0[1].scatter(Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[1]]["STEDpercent"],Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[1]]['Bleach']*100) 
seaborn.boxplot(x=Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[1]]["STEDpercent"], y=Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[1]]['Bleach']*100,ax=Ax0[1], showfliers = False, whis=(0, 100),boxprops={"facecolor": None, "edgecolor":"hotpink"},medianprops={"color": "k"})
seaborn.stripplot(x=Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[1]]["STEDpercent"], y=Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[1]]['Bleach']*100,ax=Ax0[1], size=2,color="hotpink")

seaborn.boxplot(x=Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[0]]["STEDpercent"], y=Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[0]]['Bleach']*100,ax=Ax0[0], showfliers = False, whis=(0, 100),boxprops={"facecolor": None, "edgecolor":"deepskyblue"},medianprops={"color": "k"})
seaborn.stripplot(x=Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[0]]["STEDpercent"], y=Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == labelsSPLIT[0]]['Bleach']*100,ax=Ax0[0], size=2,color="deepskyblue")




# Plot the lifetime values as a function of STED power for each sample

Ax[1,0].scatter(Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[1]]["Power"],Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[1]][ 'lifetime'])    
#Ax[1,1].scatter(Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[1]]["Power"],Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[1]]['f2'])  
Ax[1,0].scatter(Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[0]]["Power"],Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[0]][ 'lifetime'])    

#Ax[1,1].scatter(Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[0]]["Power"],Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[0]]['f2']) 

powers = numpy.unique(Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[1]]["Power"])
table_deltatau=[]
table_deltatau_median=[]
table_deltatau_std=[]
table_deltatau_sem=[]
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
    table_deltatau_median.append(numpy.median(table_deltatau_temp,axis=0))
    table_deltatau_std.append(numpy.std(table_deltatau_temp,axis=0))
    table_deltatau_sem.append(scipy.stats.sem(table_deltatau_temp,axis=0))

table_deltatau=numpy.array(table_deltatau)
table_deltatau_median=numpy.array(table_deltatau_median)
table_deltatau_std=numpy.array(table_deltatau_std)
table_deltatau_sem=numpy.array(table_deltatau_sem)
print(table_deltatau_median)
print(table_deltatau_std)
# Plot the difference in lifetime values between the 2 samples as a function of STED power
seaborn.violinplot(x=table_deltatau[:,0],y=table_deltatau[:,1],ax=Ax3,width=0.7)
tauerror=table_deltatau_std[table_deltatau_median[:,0]==0][0,1]
tau=table_deltatau_median[table_deltatau_median[:,0]==0][0,1]
Fig4, Ax4 = plt.subplots(nrows=2,ncols=3,figsize=(12,8))
Fig2, Ax2 = plt.subplots(figsize=(4,3))
Fig5, Ax5 = plt.subplots(figsize=(3,2))
Fig6, Ax6 = plt.subplots(figsize=(12,8))
Fig7, Ax7 = plt.subplots(figsize=(7,2.5))
powers=[10,20,30,40]
#powers=[5,10,15,20]

# Calculate the mean and standard deviation of the resolution metrics for the SPLIT-STED images as a function of STED power for each sample

Mean_Power_SPLIT_list=[]
STD_Power_SPLIT_list=[]
SEM_Power_SPLIT_list=[]
temp=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1])]
print(temp["lifetime"])
seaborn.boxplot(x=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1])]["Power"], y=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1])]["lifetime"],ax=Ax5, showfliers = False, whis=(0, 100),boxprops={"facecolor": None, "edgecolor":"hotpink"},medianprops={"color": "k"})
seaborn.boxplot(x=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0])]["Power"], y=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0])]["lifetime"],ax=Ax5, showfliers = False, whis=(0, 100),boxprops={"facecolor": None, "edgecolor":"deepskyblue"},medianprops={"color": "k"})

#seaborn.boxplot(Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0])], x="Power", y="lifetime",ax=Ax5)
# Add in points to show each observation
seaborn.stripplot(x=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1])]["Power"], y=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1])]["lifetime"],ax=Ax5, size=2,color="hotpink")
seaborn.stripplot(x=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0])]["Power"], y=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0])]["lifetime"],ax=Ax5, size=2,color="deepskyblue")
for i,id in enumerate(labelsSPLIT):
    Mean_Power_SPLIT=Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == id].groupby(by=["STEDpercent"]).mean(numeric_only=True).reset_index()
    STD_Power_SPLIT=Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == id].groupby(by=["STEDpercent"]).std(numeric_only=True).reset_index()
    SEM_Power_SPLIT=Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == id].groupby(by=["STEDpercent"]).sem(numeric_only=True).reset_index()
    Power_SPLIT=Overall_data_SPLIT.loc[Overall_data_SPLIT["identity"] == id].groupby(by=["STEDpercent"])

    Mean_Power_SPLIT_list.append(Mean_Power_SPLIT)
    STD_Power_SPLIT_list.append(STD_Power_SPLIT)
    SEM_Power_SPLIT_list.append(SEM_Power_SPLIT)


    for p,power in enumerate(powers):
        print(table_deltatau_std[table_deltatau_median[:,0]==power][0,1])
        xerror=table_deltatau_std[table_deltatau_median[:,0]==power][0,1]
        
        a=Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res sted stack'].item()*20

        b=STD_Power_SPLIT.loc[STD_Power_SPLIT["STEDpercent"]==power]['Res sted stack'].item()*20
        print(a,b)

        #Ax2[i,1].axhspan(a-b,a+b,alpha=0.3,facecolor=colors[i][p])
        Ax4[i,1].axhline(a,label="pSTED={}%".format(power),c=colors[i][p])


        #Ax2[i,0].errorbar(x=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],xerr=xerror,y=Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res sted stack']*20,yerr=STD_Power_SPLIT.loc[STD_Power_SPLIT["STEDpercent"]==power]['Res sted stack']*20,fmt="o",
        #           capsize=3.5, elinewidth=0.8, ms=10,ecolor="k",markeredgecolor="k",c=colors[i][p],label="pSTED={}%".format(power))
        
        Ax4[i,1].errorbar(x=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],xerr=xerror,y=Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res splitsted']*20,yerr=STD_Power_SPLIT.loc[STD_Power_SPLIT["STEDpercent"]==power]['Res splitsted']*20,fmt="o",
                   capsize=3.5, elinewidth=0.8, ms=10,ecolor="k",markeredgecolor="k",c=colors[i][p],label="pSTED={}%".format(power))
        Ax4[i,2].errorbar(x=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],xerr=xerror,y=(Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res sted stack']*20-Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res splitsted']*20)/Mean_Power_SPLIT.loc[Mean_Power_SPLIT["STEDpercent"]==power]['Res sted stack']*20,yerr=STD_Power_SPLIT.loc[STD_Power_SPLIT["STEDpercent"]==power]['Res splitsted']*20,fmt="o",
                   capsize=3.5, elinewidth=0.8, ms=10,ecolor="k",markeredgecolor="k",c=colors[i][p],label="pSTED={}%".format(power))

# Plot the confocal resolution of the 2 samples with the color of the circle indicating the lifetime difference between the 2 samples

Ax2.errorbar(x=Mean_Power_SPLIT_list[0]['Res conf stack'].mean()*20,xerr=STD_Power_SPLIT_list[0]['Res conf stack'].mean()*20,
                      y=Mean_Power_SPLIT_list[1]['Res conf stack'].mean()*20,yerr=STD_Power_SPLIT_list[1]['Res conf stack'].mean()*20,
                      fmt="o",capsize=3.5, elinewidth=0.8,ecolor="k",markeredgecolor="k",c="w",zorder=0)        
Ax2.scatter(x=Mean_Power_SPLIT_list[0]['Res conf stack'].mean()*20,
                      y=Mean_Power_SPLIT_list[1]['Res conf stack'].mean()*20,
                     c=table_deltatau_median[table_deltatau_median[:,0]==0][0,1],edgecolor="k",
                    s=100,vmin=0.75,vmax=tau,cmap="jet",label="pSTED={}%".format(0),zorder=100) 
Ax6.errorbar(x=Mean_Power_SPLIT_list[0]['Res conf stack'].mean()*20,xerr=STD_Power_SPLIT_list[0]['Res conf stack'].mean()*20,y=table_deltatau_median[table_deltatau_median[:,0]==0][0,1],yerr=table_deltatau_std[table_deltatau_median[:,0]==0][0,1],c=colors[0][0],ecolor=colors[0][0],markeredgecolor="k",fmt="o",capsize=3.5, elinewidth=0.8,zorder=0)
Ax6.errorbar(x=Mean_Power_SPLIT_list[1]['Res conf stack'].mean()*20,xerr=STD_Power_SPLIT_list[1]['Res conf stack'].mean()*20,y=table_deltatau_median[table_deltatau_median[:,0]==0][0,1],yerr=table_deltatau_std[table_deltatau_median[:,0]==0][0,1],c=colors[1][0],ecolor=colors[1][0],markeredgecolor="k",fmt="o",capsize=3.5, elinewidth=0.8,zorder=0)
Ax6.scatter(x=Mean_Power_SPLIT_list[0]['Res conf stack'].mean()*20,y=table_deltatau_median[table_deltatau_median[:,0]==0][0,1],label="pSTED={}%".format(0),c=colors[0][0],edgecolor="k",s=150)
Ax6.scatter(x=Mean_Power_SPLIT_list[1]['Res conf stack'].mean()*20,y=table_deltatau_median[table_deltatau_median[:,0]==0][0,1],label="pSTED={}%".format(0),c=colors[1][0],edgecolor="k",s=150)
    
Ax6.plot([Mean_Power_SPLIT_list[0]['Res conf stack'].mean()*20,Mean_Power_SPLIT_list[1]['Res conf stack'].mean()*20],[table_deltatau_median[table_deltatau_median[:,0]==0][0,1],table_deltatau_median[table_deltatau_median[:,0]==0][0,1]],c="k",zorder=0)
Ax7.errorbar(x=Mean_Power_SPLIT_list[0]['Res conf stack'].mean()*20,xerr=SEM_Power_SPLIT_list[0]['Res conf stack'].mean()*20,y=table_deltatau_median[table_deltatau_median[:,0]==0][0,1],yerr=table_deltatau_sem[table_deltatau_median[:,0]==0][0,1],c=colors[0][0],ecolor=colors[0][0],markeredgecolor="k",fmt="o",capsize=3.5, elinewidth=0.8,zorder=0)
Ax7.errorbar(x=Mean_Power_SPLIT_list[1]['Res conf stack'].mean()*20,xerr=SEM_Power_SPLIT_list[1]['Res conf stack'].mean()*20,y=table_deltatau_median[table_deltatau_median[:,0]==0][0,1],yerr=table_deltatau_sem[table_deltatau_median[:,0]==0][0,1],c=colors[1][0],ecolor=colors[1][0],markeredgecolor="k",fmt="o",capsize=3.5, elinewidth=0.8,zorder=0)
Ax7.scatter(x=Mean_Power_SPLIT_list[0]['Res conf stack'].mean()*20,y=table_deltatau_median[table_deltatau_median[:,0]==0][0,1],label="pSTED={}%".format(0),c=colors[0][0],edgecolor="k",s=75)
Ax7.scatter(x=Mean_Power_SPLIT_list[1]['Res conf stack'].mean()*20,y=table_deltatau_median[table_deltatau_median[:,0]==0][0,1],label="pSTED={}%".format(0),c=colors[1][0],edgecolor="k",s=75)
    
Ax7.plot([Mean_Power_SPLIT_list[0]['Res conf stack'].mean()*20,Mean_Power_SPLIT_list[1]['Res conf stack'].mean()*20],[table_deltatau_median[table_deltatau_median[:,0]==0][0,1],table_deltatau_median[table_deltatau_median[:,0]==0][0,1]],c="k",zorder=0)
  
# Plot the STED resolution of the 2 samples with the color of the circle indicating the lifetime difference between the 2 samples
for p,power in enumerate(powers):
    Ax6.errorbar(x=Mean_Power_SPLIT_list[0].loc[Mean_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,xerr=STD_Power_SPLIT_list[0].loc[STD_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,y=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],yerr=table_deltatau_std[table_deltatau_median[:,0]==power][0,1],c=colors[0][p+1],ecolor=colors[0][p+1],markeredgecolor="k",fmt="o",capsize=3.5, elinewidth=0.8,zorder=0)
    Ax6.errorbar(x=Mean_Power_SPLIT_list[1].loc[Mean_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20,xerr=STD_Power_SPLIT_list[1].loc[STD_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20,y=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],yerr=table_deltatau_std[table_deltatau_median[:,0]==power][0,1],c=colors[1][p+1],ecolor=colors[1][p+1],markeredgecolor="k",fmt="o",capsize=3.5, elinewidth=0.8,zorder=0)
    Ax6.scatter(x=Mean_Power_SPLIT_list[0].loc[Mean_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,y=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],label="pSTED={}%".format(power),c=colors[0][p+1],edgecolor="k",s=150)
    Ax6.scatter(x=Mean_Power_SPLIT_list[1].loc[Mean_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20,y=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],label="pSTED={}%".format(power),c=colors[1][p+1],edgecolor="k",s=150)
    
    Ax6.plot([Mean_Power_SPLIT_list[0].loc[Mean_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,Mean_Power_SPLIT_list[1].loc[Mean_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20],[table_deltatau_median[table_deltatau_median[:,0]==power][0,1],table_deltatau_median[table_deltatau_median[:,0]==power][0,1]],c="k",zorder=0)
    Ax7.errorbar(x=Mean_Power_SPLIT_list[0].loc[Mean_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,xerr=SEM_Power_SPLIT_list[0].loc[SEM_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,y=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],yerr=table_deltatau_sem[table_deltatau_median[:,0]==power][0,1],c=colors[0][p+1],ecolor=colors[0][p+1],markeredgecolor="k",fmt="o",capsize=3.5, elinewidth=0.8,zorder=0)
    Ax7.errorbar(x=Mean_Power_SPLIT_list[1].loc[Mean_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20,xerr=SEM_Power_SPLIT_list[1].loc[SEM_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20,y=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],yerr=table_deltatau_sem[table_deltatau_median[:,0]==power][0,1],c=colors[1][p+1],ecolor=colors[1][p+1],markeredgecolor="k",fmt="o",capsize=3.5, elinewidth=0.8,zorder=0)
    Ax7.scatter(x=Mean_Power_SPLIT_list[0].loc[Mean_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,y=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],label="pSTED={}%".format(power),c=colors[0][p+1],edgecolor="k",s=75)
    Ax7.scatter(x=Mean_Power_SPLIT_list[1].loc[Mean_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20,y=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],label="pSTED={}%".format(power),c=colors[1][p+1],edgecolor="k",s=75)
    
    Ax7.plot([Mean_Power_SPLIT_list[0].loc[Mean_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,Mean_Power_SPLIT_list[1].loc[Mean_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20],[table_deltatau_median[table_deltatau_median[:,0]==power][0,1],table_deltatau_median[table_deltatau_median[:,0]==power][0,1]],c="k",zorder=0)
    
    Ax2.errorbar(x=Mean_Power_SPLIT_list[0].loc[Mean_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,xerr=STD_Power_SPLIT_list[0].loc[STD_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,
                      y=Mean_Power_SPLIT_list[1].loc[Mean_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20,yerr=STD_Power_SPLIT_list[1].loc[STD_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20,
                      fmt="o",capsize=3.5, elinewidth=0.8,ecolor="k",markeredgecolor="k",c="w",zorder=0)        
    scat=Ax2.scatter(x=Mean_Power_SPLIT_list[0].loc[Mean_Power_SPLIT_list[0]["STEDpercent"]==power]['Res sted stack']*20,
                     y=Mean_Power_SPLIT_list[1].loc[Mean_Power_SPLIT_list[1]["STEDpercent"]==power]['Res sted stack']*20,
                     c=table_deltatau_median[table_deltatau_median[:,0]==power][0,1],edgecolor="k",
                    s=100,vmin=0.75,vmax=tau,cmap="jet",label="pSTED={}%".format(power),zorder=100) 
Fig2.colorbar(scat, ax=Ax2, location='right')
Ax2.legend()
Ax6.legend()
#Ax7.legend()
Ax7.set_xlim([80,260])
Ax7.set_ylim([0.75,0.94])
Fig0.savefig(os.path.join(savefolder,'Bleach_vs_STEDPOWER.pdf'), transparent='True', bbox_inches="tight")
Fig2.savefig(os.path.join(savefolder,'Resolution_vs_MLEmonoexp.pdf'), transparent='True', bbox_inches="tight")
Fig3.savefig(os.path.join(savefolder,'MLEmonoexp_DeltaTau.pdf'), transparent='True', bbox_inches="tight")
Fig5.savefig(os.path.join(savefolder,'MLEmonoexp_TauBoxplot.pdf'), transparent='True', bbox_inches="tight")
Fig6.savefig(os.path.join(savefolder,'MLEmonoexp_DeltaTau_simple.pdf'), transparent='True', bbox_inches="tight")
Fig7.savefig(os.path.join(savefolder,'MLEmonoexp_DeltaTau_simple_sem.pdf'), transparent='True', bbox_inches="tight")
# Plot the lifetime values as a function of STED power for each sample
seaborn.boxplot(x=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1])]["Power"], y=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1])]["lifetime"],ax=Ax5, showfliers = False, whis=(0, 100),boxprops={"facecolor": None, "edgecolor":"hotpink"},medianprops={"color": "k"})
seaborn.boxplot(x=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0])]["Power"], y=Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0])]["lifetime"],ax=Ax5, showfliers = False, whis=(0, 100),boxprops={"facecolor": None, "edgecolor":"deepskyblue"},medianprops={"color": "k"})

# Statistical tests to compare the lifetime values of the 2 samples for each STED power
print("################################################################")
print("STATS Tau PSD95")

stat, p = scipy.stats.shapiro(Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1])]["lifetime"] )
print("Shapiro result",stat, p)
if p<0.05:
    print("All together,PSD95 lifetimes are not normal")

stat, p = scipy.stats.shapiro(Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0])]["lifetime"] )
print("Shapiro result",stat, p)
if p<0.05:
    print("All together,Bassoon lifetimes are not normal")

powers = numpy.unique(Overall_data_MLE.loc[Overall_data_MLE["identity"] == labelsSPLIT[1]]["Power"])
for power in powers:
    stat, p = scipy.stats.shapiro(Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1])&(Overall_data_MLE["Power"]==power)]["lifetime"])
    print("Shapiro result PSD95",stat, p)
    if p<0.05:
        print("Taus are not normal")
    stat, p = scipy.stats.shapiro(Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0])&(Overall_data_MLE["Power"]==power)]["lifetime"])
    print("Shapiro result Basoon",stat, p)
    if p<0.05:
        print("Taus are not normal")
try:
    stat, p =scipy.stats.kruskal(*[Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1])&(Overall_data_MLE["Power"]==power)]["lifetime"] for power in powers])
    print("kruskal result PSD95",stat, p)
    stat, p =scipy.stats.kruskal(*[Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0])&(Overall_data_MLE["Power"]==power)]["lifetime"] for power in powers])
    print("kruskal result Bassoon",stat, p)
except ValueError:
    print("Kruskal did not work")
    
print("get_Significance PSD95")
sig=get_significance([Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[1])&(Overall_data_MLE["Power"]==power)]["lifetime"] for power in powers], verbose=True)
print(sig)
print("get_Significance Bassoon")
sig=get_significance([Overall_data_MLE.loc[(Overall_data_MLE["identity"] == labelsSPLIT[0])&(Overall_data_MLE["Power"]==power)]["lifetime"] for power in powers], verbose=True)
print(sig)
print("###############################################################")
plt.show()