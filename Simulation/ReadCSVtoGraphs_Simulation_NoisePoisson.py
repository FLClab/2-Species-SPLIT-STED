"""
This script generates performance graphs based on simulation results with added synthetic background. It reads CSV files generated
 by Simulation_2SpeciesSPLIT-STED_looppowers_NoiseonMixed_Poisson.py containing metrics of unmixed simulation images, such as resolution and SQUIRREL.
 The script groups the metrics by condition such as depletion power and noise type and plots the results using matplotlib and seaborn. 
 
    It also generates histograms of the Poisson noise distributions used in the simulations.
"""

import os
import glob

import numpy
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn
import easygui
import itertools
from sys import path as path1; 
Functionspath=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Functions")
path1.append(Functionspath)
from statistics_functions import get_significance
import scipy

matplotlib.rcParams['axes.linewidth'] = 0.8

fig,ax=plt.subplots(figsize=(3,1))
lambdas=[3,5,7]
for l in lambdas:
    noisegrid=numpy.random.poisson(lam=l, size=10000).astype(int)
    ax.hist(noisegrid, bins=21, range=(0, 20),histtype='step',density=True, fill=False, label=f'λ={l}')
ax.legend()
#plt.show()    
ax.set_xlabel('Photon Count')
ax.set_ylabel('Frequency')



# Ask user name of the folder where the plots will be saved and create it
savefoldername =str(input("Name of folder: "))
savefolder=os.path.join(os.path.expanduser("~/Desktop"),savefoldername + "_Simulation_Noise")
os.makedirs(savefolder,exist_ok=True)

#Path to folder containing the .csv files
#csvpath = easygui.fileopenbox(default=os.path.expanduser("~/Desktop"),multiple=False)
folder=easygui.diropenbox(default=os.path.expanduser("~/Desktop"))
folders=[folder]



colors=[["#7ce8ffff","#55d0ffff","#00acdfff","#2895cbff","#1676a9ff","#00456bff"],["#fcbcd7ff","#f9a3cbff","#ef87beff","#e569b3ff","#bf4290ff","#5d1a43ff","#981164ff"]]

# Read the .csv files and concatenate them
cumdf=[]
csvfull=[]
# Loop over the folders containing the CSV files and pool the data into a single dataframe.
for i,folder in enumerate(folders):

    csvlist=glob.glob(os.path.join(folder,"*.csv"))
    print(csvlist)


    Overall_data= pd.concat(map(pd.read_csv, csvlist))
    cumdf.append(Overall_data)

    print(Overall_data.shape)
    print(list(Overall_data.columns))
    print(Overall_data)

Overall_data.dropna(axis=0, how='all', inplace=True)
for path in Overall_data["image1"]:
  
    Overall_data["image1"]= os.path.basename(path)
for path in Overall_data["image2"]:
    Overall_data["image2"]= os.path.basename(path)

#Overall_data= pd.read_csv(csvpath)
#print(Overall_data.shape)
#print(list(Overall_data.columns))

df=Overall_data[Overall_data["Noise type"]=="Uniform"]
df=df[df["Noise number photons"]==0]

dfcontrol=df.copy()
powersC=numpy.unique(df["Power"])
print(powersC)
SQ1C=numpy.array(df[["squirrel_f1"]])
SQ2C=numpy.array(df[["squirrel_f2"]])
PowC=numpy.array(df[["Power"]])

dfmean=df.groupby("Power")[ ['Power', 'squirrel_f1', 'squirrel_f2']].mean(numeric_only=True)
dfstd=df.groupby("Power")[ ['Power',  'squirrel_f1', 'squirrel_f2']].sem(numeric_only=True)

df=Overall_data[Overall_data["Noise type"]=="2Species"]
df2s=df.copy()
dfmean2s=df.groupby("Power")[ ['Power', 'squirrel_f1', 'squirrel_f2']].mean(numeric_only=True)
dfstd2s=df.groupby("Power")[ ['Power',  'squirrel_f1', 'squirrel_f2']].sem(numeric_only=True)
mergedf=pd.merge(dfcontrol,df2s, on=["Power","image1","image2"], suffixes=('_control', '_noise'))
mergedf["Delta_SQ1"]=mergedf["squirrel_f1_noise"]-mergedf["squirrel_f1_control"]
mergedf["Delta_SQ2"]=mergedf["squirrel_f2_noise"]-mergedf["squirrel_f2_control"]
mergedfmean2s=mergedf.groupby("Power")[ ["Power","Delta_SQ1", "Delta_SQ2"]].mean(numeric_only=True)
mergedfstd2s=mergedf.groupby("Power")[ ["Power","Delta_SQ1", "Delta_SQ2"]].sem(numeric_only=True)

types=["Uniform","IRF","Alexa647"]
fig2,ax2=plt.subplots(ncols=3,nrows=2,sharex=True,sharey=True,figsize=(10,6))
for t,type in enumerate(types):
    fig3,ax3=plt.subplots(ncols=1,nrows=2,sharex=True,sharey=True,figsize=(4,8))
    #fig1,ax1=plt.subplots(figsize=(8,6))
    df=Overall_data[Overall_data["Noise type"]==type]
    pixels=numpy.unique(df["Noise percent pixels"])
    photons=numpy.unique(df["Noise number photons"])
    print(type)
    print(pixels)
    print(photons)
    if 0 in pixels and t<3:
        print("Removing 0 from pixels")
        pixels=numpy.delete(pixels,numpy.where(pixels==0)[0][0])
    print(pixels)
    if 0 in photons and t<3:
        print("Removing 0 from photons")
        photons=numpy.delete(photons,numpy.where(photons==0)[0][0])
    print(photons)
   
    for i,p in enumerate(pixels):
        if p==0 and t<3:
            continue
        df2=df[df["Noise percent pixels"]==p]

        #ax2[t,i].errorbar(powersC,dfmean["squirrel_f1"],yerr=dfstd["squirrel_f1"],c="k",fmt="o", label="Control F1",ecolor="k", capsize=5, elinewidth=0.8, ms=5)
        #ax2[t,i].errorbar(powersC,dfmean["squirrel_f2"],yerr=dfstd["squirrel_f2"],c="grey", fmt="o", label="Control F2",ecolor="grey", capsize=5, elinewidth=0.8, ms=5)
        ax3[0].errorbar(powersC,dfmean["squirrel_f1"],yerr=dfstd["squirrel_f1"],c="k",fmt="o", label="Control F1",ecolor="k", capsize=5, elinewidth=0.8, ms=5)
        ax3[1].errorbar(powersC,dfmean["squirrel_f2"],yerr=dfstd["squirrel_f2"],c="grey", fmt="o", label="Control F2",ecolor="grey", capsize=5, elinewidth=0.8, ms=5)

        ax2[0,t].fill_between(powersC, dfmean["squirrel_f1"]-dfstd["squirrel_f1"], dfmean["squirrel_f1"]+dfstd["squirrel_f1"], color="k", alpha=0.2)
        ax2[1,t].fill_between(powersC, dfmean["squirrel_f2"]-dfstd["squirrel_f2"], dfmean["squirrel_f2"]+dfstd["squirrel_f2"], color="k", alpha=0.2)
        ax2[0,t].plot(powersC,dfmean["squirrel_f1"],c="k", label="2 Species SPLIT-STED \n No Background - F1", marker="o", ls='--')
        ax2[1,t].plot(powersC,dfmean["squirrel_f2"],c="k", label="2 Species SPLIT-STED \n No Background - F2", marker="o", ls='--')

        ax2[0,t].fill_between(powersC, dfmean2s["squirrel_f1"]-dfstd2s["squirrel_f1"], dfmean2s["squirrel_f1"]+dfstd2s["squirrel_f1"], color="grey", alpha=0.2)
        ax2[1,t].fill_between(powersC, dfmean2s["squirrel_f2"]-dfstd2s["squirrel_f2"], dfmean2s["squirrel_f2"]+dfstd2s["squirrel_f2"], color="grey", alpha=0.2)
        ax2[0,t].plot(powersC,dfmean2s["squirrel_f1"],c="grey", label="2 Species STED-FLIM\n No Background - F1", marker="*", ls='--')
        ax2[1,t].plot(powersC,dfmean2s["squirrel_f2"],c="grey", label="2 Species STED-FLIM \n No Background - F2", marker="*", ls='--')

        #ax2[0,t].fill_between(powersC, mergedfmean2s["Delta_SQ1"]-mergedfstd2s["Delta_SQ1"], mergedfmean2s["Delta_SQ1"]+mergedfstd2s["Delta_SQ1"], color="k", alpha=0.2)
        #ax2[1,t].fill_between(powersC, mergedfmean2s["Delta_SQ2"]-mergedfstd2s["Delta_SQ2"], mergedfmean2s["Delta_SQ2"]+mergedfstd2s["Delta_SQ2"], color="k", alpha=0.2)
        #ax2[0,t].plot(powersC,mergedfmean2s["Delta_SQ1"],c="k", label="2 Species STED-FLIM \n No Background - F1", marker="o", ls='--')
        #ax2[1,t].plot(powersC,mergedfmean2s["Delta_SQ2"],c="k", label="2 Species STED-FLIM \n No Background - F2", marker="o", ls='--')

        for j,n in enumerate(photons):

            df3=df2[df2["Noise number photons"]==n]
            powers=numpy.unique(df3["Power"])
            print(powers)
            #powers=[10,20,30]
            mergedf=pd.merge(dfcontrol,df3, on=["Power","image1","image2"], suffixes=('_control', '_noise'))
            mergedf["Delta_SQ1"]=mergedf["squirrel_f1_noise"]-mergedf["squirrel_f1_control"]
            mergedf["Delta_SQ2"]=mergedf["squirrel_f2_noise"]-mergedf["squirrel_f2_control"]
            mergedfmean=mergedf.groupby("Power")[ ["Power","Delta_SQ1", "Delta_SQ2"]].mean(numeric_only=True)
            mergedfstd=mergedf.groupby("Power")[ ["Power","Delta_SQ1", "Delta_SQ2"]].sem(numeric_only=True)
            print(mergedf)


            df3mean=df3.groupby("Power")[ ['Power', 'squirrel_f1', 'squirrel_f2']].mean(numeric_only=True)
            df3std=df3.groupby("Power")[ ['Power', 'squirrel_f1', 'squirrel_f2']].sem(numeric_only=True)

      

        

            #ax2[t,i].scatter(df3["Power"],df3["squirrel_f1"],c=colors[0], label="{} photons".format(n))
            #ax2[t,i].scatter(df3["Power"],df3["squirrel_f2"],c=colors[1], label="{} photons".format(n))
            #ax2[t,i].axhline(0,c="k",ls="--")
            #ax2[0,t].scatter(mergedf["Power"],mergedf["Delta_SQ1"],c=colors[0][j], label="λ={} F1".format(n), alpha=0.5, s=10)
            #ax2[1,t].scatter(mergedf["Power"],mergedf["Delta_SQ2"],c=colors[1][j], label="λ={} F2".format(n), alpha=0.5, s=10)

            #ax2[0,t].fill_between(powers, mergedfmean["Delta_SQ1"]-mergedfstd["Delta_SQ1"], mergedfmean["Delta_SQ1"]+mergedfstd["Delta_SQ1"], color=colors[0][j], alpha=0.2)
            #ax2[1,t].fill_between(powers, mergedfmean["Delta_SQ2"]-mergedfstd["Delta_SQ2"], mergedfmean["Delta_SQ2"]+mergedfstd["Delta_SQ2"], color=colors[1][j], alpha=0.2)
            #ax2[0,t].plot(powers,mergedfmean["Delta_SQ1"],c=colors[0][j], label="λ={} F1".format(n), marker="o", ls='--')
            #ax2[1,t].plot(powers,mergedfmean["Delta_SQ2"],c=colors[1][j], label="λ={} F2".format(n), marker="o", ls='--')

            ax2[0,t].fill_between(powers, df3mean["squirrel_f1"]-df3std["squirrel_f1"], df3mean["squirrel_f1"]+df3std["squirrel_f1"], color=colors[0][j], alpha=0.2)
            ax2[1,t].fill_between(powers, df3mean["squirrel_f2"]-df3std["squirrel_f2"], df3mean["squirrel_f2"]+df3std["squirrel_f2"], color=colors[1][j], alpha=0.2)
            ax2[0,t].plot(powers,df3mean["squirrel_f1"],c=colors[0][j], label="λ={} F1".format(n), marker="o", ls='--')
            ax2[1,t].plot(powers,df3mean["squirrel_f2"],c=colors[1][j], label="λ={} F2".format(n), marker="o", ls='--')


            #ax2[0,t].errorbar(powers,mergedfmean["Delta_SQ1"],yerr=mergedfstd["Delta_SQ1"],fmt="o",c=colors[0][j], label="λ={} F1".format(n),ecolor=colors[0][j], capsize=10, elinewidth=0.8, ms=8)
            #ax2[1,t].errorbar(powers,mergedfmean["Delta_SQ2"],yerr=mergedfstd["Delta_SQ2"],fmt="o",c=colors[1][j], label="λ={} F2".format(n),ecolor=colors[1][j], capsize=10, elinewidth=0.8, ms=8)
            #ax2[t,i].errorbar(powers,df3mean["squirrel_f1"],yerr=df3std["squirrel_f1"],fmt="o",c=colors[0][j], label="{} photons F1".format(n),ecolor=colors[0][j], capsize=5, elinewidth=0.8, ms=5)
            #ax2[t,i].errorbar(powers,df3mean["squirrel_f2"],yerr=df3std["squirrel_f2"],fmt="o",c=colors[1][j], label="{} photons F2".format(n),ecolor=colors[1][j], capsize=5, elinewidth=0.8, ms=5)
            ax3[0].errorbar(powers,df3mean["squirrel_f1"],yerr=df3std["squirrel_f1"],fmt="o",c=colors[0][j], label="λ={} F1".format(n),ecolor=colors[0][j], capsize=5, elinewidth=0.8, ms=5)
            ax3[1].errorbar(powers,df3mean["squirrel_f2"],yerr=df3std["squirrel_f2"],fmt="o",c=colors[1][j], label="λ={} F2".format(n),ecolor=colors[1][j], capsize=5, elinewidth=0.8, ms=5)


        ax2[0,t].set_title("{}".format(type))
        ax2[0,t].set_xlim([5,45])
        ax2[1,t].set_xlim([5,45])
        #if t<2:
        #    ax2[t,i].set_xticklabels([])
        #if i>0:
        #    ax2[t,i].set_yticklabels([])
        ax3[i].set_title("{}".format(type))
        ax2[1,t].set_xlabel("Power")
        ax2[0,0].set_ylabel("SQUIRREL F1")
        ax2[1,0].set_ylabel("SQUIRREL F2")
        ax3[0].set_xlabel("Power")
        ax3[0].set_ylabel("SQUIRREL F1")
        ax3[1].set_xlabel("Power")
        ax3[1].set_ylabel("SQUIRREL F2")
        ax3[0].legend(frameon=False)
        ax3[1].legend(frameon=False)
        ax2[0,1].legend(frameon=False,loc='upper right')
        ax2[1,1].legend(frameon=False,loc='upper right')

    fig3.savefig(os.path.join(savefolder,"SQUIRREL_{}.pdf".format(type)),transparent=True, bbox_inches='tight')
fig2.subplots_adjust(wspace=0,hspace=0)
fig2.savefig(os.path.join(savefolder,"DeltaSQUIRREL.pdf"),transparent=True, bbox_inches='tight')
fig.savefig(os.path.join(savefolder,"Poisson_Dists.pdf"),transparent=True, bbox_inches='tight')
plt.show()
