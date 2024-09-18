import os
import functools
import warnings
import time
import csv
import matplotlib; matplotlib.use("TkAgg")
from matplotlib import pyplot
import sys
sys.path.append(os.path.join("C:", os.sep, "Users", "abberior","Documents","GitHub","Abberior-STED"))
import skimage.io
# from skimage.external.tifffile import TiffWriter
import yaml
from abberior import microscope, user, utils
from datetime import date
import pandas
import numpy
import random

#liste = numpy.array(pandas.read_csv(r"C:\Users\abberior\Documents\GitHub\RESOLFT-Abberior\rsfp_rsGEFP2.csv")).transpose()[2:,:]
#powers=numpy.linspace(5,40,num=5,dtype=int).tolist()
powers = [5,10,15,20]
powers = [10,20,30,40]
powers = list(numpy.repeat(powers, 5))

random.shuffle(powers)
print(powers)
# power = [10,20,30,40,50]

n = len(powers)
overview='561' #À changer avec les rsfps rouges
today = date.today()
folder = os.path.join("C:", os.sep, "Users", "abberior", "Desktop","DATA","Andreanne",
                      "{}-{:02d}-{:02d}-FLIM_Power_Bleach_DifferentRegions_".format(today.year, today.month, today.day))
id = input('Enter sample region identifier: ')
output = os.path.join(folder, id)

#imsavepath = os.path.join(output, "Images")
configsavepath = os.path.join(output, "Configurations")

#os.makedirs(imsavepath, exist_ok=True)
os.makedirs(configsavepath, exist_ok=True)
config_overview = microscope.get_config("Setting overview configuration.")
config_Pre = microscope.get_config("Setting Confocal before configuration.")
config_Read = microscope.get_config("Setting STED configuration.")
#config_Read2 = microscope.get_config("Setting HighPower STED configuration.")
config_Post = microscope.get_config("Setting Confocal after configuration.")
config_window=config_Read
configs=[config_Pre,config_Read,config_Post]
dictionnary = {}


t=0
while n > 0:
    #t = int(input("Vous êtes rendus à quelle région?")) #À CHANGER SELON LE NUMÉRO DE LA MESURE À FAIRE, CAR ON VEUT DES NOMS DE RÉGIONS DIFFÉRENTS
    iiii = input('Find a region then press enter')
    if iiii == 'non':
        break
    rect_region_offset, regions, rectangles = user.get_rect_regions(config_overview=config_overview,overview=overview)

    # For each rectangle selected by user, set imaging parameters, acquire image and save configuration and images
    for (x, y), (width, height) in zip(rect_region_offset, regions):

        # This will output a .csv file for every region that were imaged
        # parameters x,y are the imaging region offset center
        # width, height are the imaging size
        with open(os.path.join(configsavepath, "imaging_parameters_Region{}".format(t)), "a") as csvfile:
            writer = csv.writer(csvfile, delimiter="\t")
            writer.writerow([x, y,width,height])

        print('regions saved to csv')

        for config_window in configs:
            # Set microscope configurations to image selected region
            microscope.set_imagesize(config_window, width, height)
            microscope.set_offsets(config_window, x, y)

        microscope.set_power(config_Read, powers[0], 6, channel_id=0) #id pour le 775
        print('microscope parameters set')


        stacks, _ = microscope.acquire_multi_saveasmsr(configs, os.path.join(output, id + '_{}_'.format(t)+str(powers[0])+'percentSTED.msr'))

        del powers[0]

        t+=1
        if powers == []:
            print('ALL DONE')
            n = len(powers)
            break

    n = len(powers)
    print(f"Voici le nombre de régions prises à date: {t}")

    with open(os.path.join(configsavepath, "imspector_config_window_{}".format(t)), "w") as f:
        config = config_Read.parameters("")
        yaml.dump(config, f)
