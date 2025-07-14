
"""
This script automates the acquisition of STED-FLIM images with varying depletion powers using an Abberior STED microscope.
The script repeats and shuffles a list of depletion powers to be used for the acquisition to ensure randomization.
It iterates through the list of powers and asks the user to select regions to be imaged. For each region, the script sets
imaging parameters,acquires images using the established sequence, and saves the configuration and images.
"""

import os
import csv
import matplotlib; matplotlib.use("TkAgg")
from matplotlib import pyplot
import yaml
from abberior import microscope, user, utils
from datetime import date
import numpy
import random

# List of powers to be used for the acquisition, in % and the number of repetitions
#powers = [5,10,15,20]
powers = [10,20,30,40]
reps = 5
# Repeat the power list and shuffle the list
powers = list(numpy.repeat(powers, reps))
random.shuffle(powers)
print(powers)
n = len(powers)

# Create folder to save images and configurations
overview='561'
today = date.today()
folder = os.path.join("C:", os.sep, "Users", "abberior", "Desktop","DATA","Andreanne",
                      "{}-{:02d}-{:02d}-FLIM_Power_Bleach_DifferentRegions_".format(today.year, today.month, today.day))
id = input('Enter sample region identifier: ')
output = os.path.join(folder, id)
configsavepath = os.path.join(output, "Configurations")
os.makedirs(configsavepath, exist_ok=True)

# Ask the user to select the measurement window corresponding to each imaging step
config_overview = microscope.get_config("Setting overview configuration.")
config_Pre = microscope.get_config("Setting Confocal before configuration.")
config_Read = microscope.get_config("Setting STED configuration.")
config_Post = microscope.get_config("Setting Confocal after configuration.")
config_window=config_Read

configs=[config_Pre,config_Read,config_Post]# List of measurements to be run in order


t=0
while n > 0:

    iiii = input('Find a region then press enter')
    if iiii == 'non':
        break
    rect_region_offset, regions, rectangles = user.get_rect_regions(config_overview=config_overview,overview=overview)

    # For each rectangle selected by user, set imaging parameters, acquire image and save configuration and images
    for (x, y), (width, height) in zip(rect_region_offset, regions):

        # This will output a .csv file for every region that is imaged
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
        # Set the depletion power to the first value in the list
        microscope.set_power(config_Read, powers[0], 6, channel_id=0) #id for the 775 laser
        print('microscope parameters set')

        # Acquire images and save them
        stacks, _ = microscope.acquire_multi_saveasmsr(configs, os.path.join(output, id + '_{}_'.format(t)+str(powers[0])+'percentSTED.msr'))

        del powers[0]

        t+=1
        if powers == []:
            print('ALL DONE')
            n = len(powers)
            break

    n = len(powers)
    print(f"Here is the number of regions that have been imaged: {t}")

# Save the configuration used for the acquisition to a YAML file
    with open(os.path.join(configsavepath, "imspector_config_window_{}".format(t)), "w") as f:
        config = config_Read.parameters("")
        yaml.dump(config, f)
