""" Computes IRF (FWHM) value with a gaussian fit
    Use a FLIM confocal gold bead images because the 
    gold reflects all laser beam in the detector
"""

import matplotlib.pyplot as plt
import glob
import numpy
from matplotlib_scalebar.scalebar import ScaleBar
import easygui

from scipy.optimize import curve_fit
from scipy import asarray as ar,exp


import os.path
from sys import path as path1; path1.append('/Users/marielafontaine/Documents/GitHub/Abberior-STED-FLIM/Functions')
dossier = os.path.expanduser("~/Documents/Github/Abberior-STED-FLIM/Functions")
path1.append(dossier)
from Main_functions import (choose_msr_file, get_foreground)
# -----------------------------------------------------------
#                  Sélectionner les images
#easygui.diropenbox(default=os.path.expanduser("~Desktop"))

filename = '/Users/marielafontaine/valeria-s3/flclab-abberior-sted/mlafontaine/22-06-23_Gold_Bead'
mapcomp = {'CONF' : 'STAR 635P_CONF {0}'}



path=os.path.join(filename,"*.msr")
images = glob.glob(path)
print('There are ',len(images), 'Images in this folder')


# -----------------------------------------------------------
#     Choisir un fichier msr parmi plusieurs

imagemsr = choose_msr_file(images)
print(imagemsr.keys())


#plt.style.use('dark_background')

for key in mapcomp :
    image1 = imagemsr[mapcomp[key]]

    dim = image1.shape
    if dim[2] > 250 :
        image1 = image1[:,:,0:250]
        dim=image1.shape

    imsum= numpy.sum(image1, axis=2)
    fig, ax = plt.subplots()
    ax.imshow(imsum)
    ax.axis('off')
    scalebar = ScaleBar(0.02, "um", length_fraction=0.25)
    ax.add_artist(scalebar)
    imgplot1 = ax.imshow(imsum, cmap='hot')
    cbar =fig.colorbar(imgplot1)

    seuil = get_foreground(imsum)

    fig_seuil, ax_seuil = plt.subplots()
    ax_seuil.imshow(imsum > seuil) 
    scalebar = ScaleBar(0.02, "um", length_fraction=0.25)
    #ax.add_artist(scalebar)
    ax_seuil.axis('off')


    # -----------------------------------------------------------
    fig_gaus, ax_gaus = plt.subplots()

    y = image1[imsum > seuil, :] 

    y= numpy.sum(y, axis=0)
    y1 = y[ : ] 
    hist1 = y / y.max()
    y= y1 / y1.sum()

    absci = numpy.linspace(0,y.shape[0]-1, num =y.shape[0])*0.08
    hist = ax_gaus.bar(absci, hist1, width=0.08)

    def gaus(x,a,x0,sigma):
        return a*numpy.exp(-(x-x0)**2/(2*sigma**2)) 

    popt,pcov = curve_fit(gaus,absci,y)
    gaussian = gaus(absci,*popt) / gaus(absci,*popt).max()

    sigma = popt[2]
    mu = popt[1]

    FWHM = 2.3548 * sigma
    print('FWHM', FWHM)
    print('sigma', sigma)
    print('mu', mu)

    y /= y.max()
    ax_gaus.plot(absci, y , 'b+:', label='DATA')
    ax_gaus.plot(absci,gaussian,'ro:',label='fit')
    ax_gaus.set_xlabel('time [ns]', fontsize = 20)
    ax_gaus.set_ylabel('Normalised intensity', fontsize = 20)
    ax_gaus.xaxis.set_tick_params(labelsize=20)
    ax_gaus.yaxis.set_tick_params(labelsize=20)

    ax_gaus.legend(fontsize = 20)
    plt.show()

