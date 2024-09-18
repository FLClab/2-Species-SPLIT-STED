""" Prend un dossier qui contient plusieurs fichiers msr
    Ouvre chaque fichier et vérifie si les images indiquées par keys ont plus de photons que les valeurs dans limits
"""


import os
import glob
import numpy
import matplotlib.pyplot as plt
import easygui 
from sys import path as syspath; 
dossier = os.path.expanduser("~/Documents/Github/2-Species-SPLIT-STED/Functions")
syspath.append(dossier)
from Main_functions import load_msr 

filename = easygui.diropenbox(default=os.path.expanduser("~Desktop"))

#keys=["STED_635P {2}","Conf_635P {2}"]
keys=['STED 561 {11}','Confocal_561 {11}']
#keys=['STED_594 {2}','Conf_ 594 {2}']
#keys=['STED640 {10}' ,'Conf640 {10}']#Live Cy5
############Limits are chosen to respect a max 1% probability of 2 photons/pulse (5 photons/us) #######################################
#limits=[450,225] #Fixed Cy5 with and without 60P control
#limits=[450,75] #Fixed Cy5 and Cy3 24/08
#limits=[560,210] #Fixed Cy3 with 50P control
#limits=[600,300] #Fixed Cy3 without 50P control
limits=[225,150] #Live Cy3
limits=[456,120]
#limits=[1200,225] #Live Cy5
#filename = '/Users/marielafontaine/valeria-s3/flclab-abberior-sted/mlafontaine/2022-08-24-FLIM_Power_Bleach/Tetraspeck_5-to-40percentSTED_2conf-1'


path=os.path.join(filename,"*.msr")
images = glob.glob(path)
print('There are ',len(images), 'Images in this folder')


for i, imagei in enumerate(images):
    
    imagemsr = load_msr(imagei)

    #print(imagemsr.keys())
    for k,key in enumerate(keys) :
        #value = mapcomp[key]
        image1=imagemsr[key]
        #print(i, image1.shape, key)
        print(os.path.basename(imagei)) 
        if len(image1.shape) == 3 :
            imsum= numpy.sum(image1, axis=2)
            print(key,numpy.max(imsum))
            if numpy.max(imsum)>limits[k]:
                print("")
                #print(key, "  We should consider removing this image it has: ",numpy.count_nonzero(imsum>limits[k]), 'pixels above the limit with a max of',numpy.max(imsum))
        elif len(image1.shape) == 1 :
            continue

    
            



# cbar =fig.colorbar(imgplot1)
# scalebar = ScaleBar(0.02, "um", length_fraction=0.25)
# ax.add_artist(scalebar)
#cbar.set_label("Intensité")


