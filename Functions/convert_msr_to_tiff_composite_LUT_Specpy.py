# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import skimage
from specpy import File
from skimage import io
import os
import glob
import numpy
import tifffile
from tiffwrapper import imwrite
import easygui

def load_msr(msrPATH):
    """Loads the msr data from the msr file. The data is in a numpy array.

    :param msrPATH: The path to the msr file

    :returns: A dict with the name of the stack and the corresponding data
    """
    imageObj = File(msrPATH, File.Read)
    outputDict = {}
    outputDictMeta = {}
    for stack in range(imageObj.number_of_stacks()):
        stackObj = imageObj.read(stack) # stack object
        metadata=stackObj.parameters("doc")
        metadata=metadata['ExpControl']['scan'][ 'range']
        outputDict[stackObj.name()] = stackObj.data().squeeze() # removes shape of 1
        outputDictMeta[stackObj.name()]=metadata
    return outputDict,outputDictMeta

def load_msr_nometa(msrPATH):
    """Loads the msr data from the msr file. The data is in a numpy array.

    :param msrPATH: The path to the msr file

    :returns: A dict with the name of the stack and the corresponding data
    """
    imageObj = File(msrPATH, File.Read)
    outputDict = {}

    for stack in range(imageObj.number_of_stacks()):
        stackObj = imageObj.read(stack) # stack object

        outputDict[stackObj.name()] = stackObj.data().squeeze() # removes shape of 1

    return outputDict

##When the code is run, a navigator window will appear. Browse to and select folder containing .msr images to convert
print("Press Alt+tab to find the windows browser \n")
filename = easygui.diropenbox(default=os.path.expanduser("~Desktop"))
print(filename)

path=os.path.join(filename,"*.msr")

images = glob.glob(path)

print('There are ',len(images), 'Images in this folder')

# Image identifiers (Channel names) to be included in images. List multiple channels when a composite is wanted
mapcomp={
         '': ["Confocal_561 {11}", "STED 561 {11}"]
         }
mapcomp={ 
        '': ['Conf_635P {2}', 'STED_635P {2}']
        }

mapfoldernames={}
for key in mapcomp:
    outpath=os.path.join(os.path.expanduser("~/Desktop"),os.path.basename(filename))
    os.makedirs(outpath, exist_ok=True)
    mapfoldernames[key]=outpath

for c,imagei in enumerate(images):
    print("{}/{}".format(c,len(images)),imagei)
    try:
        meta=True
        imagemsr,metadata=load_msr(imagei)

    except RuntimeError:
        try:
            imagemsr=load_msr_nometa(imagei)
            meta=False
        except RuntimeError:
            print('This file is busy, try again later')
            continue
    for key in mapcomp:
        outpath=os.path.join(filename,key)

        channels=mapcomp[key]

        if len(channels)==1:
            try:
                image1=imagemsr[channels[0]]
                #print(metadata[channels[0]])
                if meta:
                    pixelX=metadata[channels[0]]['x']['psz']
                    pixelY = metadata[channels[0]]['y']['psz']
                image1=numpy.expand_dims(image1, axis=0)
                image1=numpy.moveaxis(image1,3,0)
                print(image1.shape)
                filenameout=os.path.join(mapfoldernames[key],os.path.basename(imagei).split(".msr")[0]+"{}.tiff".format(key))
                if meta:
                    imwrite(file=filenameout, data=image1.astype(numpy.uint16), composite=True,pixelsize=(pixelX * 1e+6,pixelY* 1e+6))
                else:
                    imwrite(file=filenameout, data=image1.astype(numpy.uint16), composite=True)


            except KeyError:
                print('No such key', channels[0])
                print(imagemsr.keys())
                continue
            except RuntimeError:
                print('Corrupted file, skipping...')
                continue


            #io.imsave(imagei.split(".msr")[0]+"_{}.tiff".format(key),image1)

        elif len(channels)>1:
            try:
                image1=[imagemsr[channel] for channel in channels]

                if meta:
                    pixelX=metadata[channels[0]]['x']['psz']
                    pixelY = metadata[channels[0]]['y']['psz']


            except KeyError:
                print('No such key', channels)
                print(imagemsr.keys())
                continue
            except RuntimeError:
                print('Corrupted file, skipping...')
                continue

            imagecomp=numpy.stack(image1)
            #print(imagecomp.shape)
            imagecomp=numpy.moveaxis(imagecomp,3,0)
            print(imagecomp.shape)


            filenameout=os.path.join(mapfoldernames[key],os.path.basename(imagei).split(".msr")[0]+"{}.tiff".format(key))
            if meta:
                imwrite(file=filenameout, data=imagecomp.astype(numpy.uint16), composite=True,pixelsize=(pixelX * 1e+6,pixelY* 1e+6))
            else:
                imwrite(file=filenameout, data=imagecomp.astype(numpy.uint16), composite=True)
