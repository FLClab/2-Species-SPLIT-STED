
from pylab import *
from statsmodels.tsa.stattools import acf
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib
matplotlib.use('TkAgg')
import numpy
import numpy as np

from matplotlib import gridspec, pyplot
from skimage import measure,io
#import beadlines

import os
import glob
import csv
import tifffile
from tiffwrapper import make_composite
import easygui
from specpy import File
#from convertmsr_bioformatsAB import MSRReader

rows= (0,1,2,0,1,2)
columns= (0,0,0,1,1,1)
plotcolors = ( 'gray','blue','gold','gray','blue','gold')
labels=["2 species phasor-FLIM - Confocal","2 species phasor-FLIM - STED","2 species SPLIT-STED","2 species phasor-FLIM - Confocal","2 species phasor-FLIM - STED","2 species SPLIT-STED"]
luts=("Cyan Hot", "Magenta Hot")
rangelut=((0, 100), (0, 100))
dispchannel=[2,5]
Width=5
pixelsize=20


class BeadlinesFigure:
    """Class that displays an image in a matplotlib window and allows the user to draw line segments.
    Line coordinates are saved in a list.
    Takes input image data as numpy array as parameter."""

    def __init__(self, img_data):

        self.fig, self.ax = plt.subplots(figsize=(15, 8))
        self.clicked_once = False
        self.xs = []
        self.ys = []
        self.linelist = []
        self.current_marker = None
        self.img_data = img_data

        plt.imshow(img_data, cmap='hot')
        #plt.colorbar()

        self.ax.set_title('Draw lines by clicking on endpoints. r removes last line ')

        self.fig.canvas.mpl_connect('button_press_event', self.click_event)
        self.fig.canvas.mpl_connect('key_press_event', self.remove_selection)
        self.fig.canvas.mpl_connect('close_event', self.get_values)

        plt.show()

    def click_event(self, event):
        """Allows the drawing of lines on the figure. First click saves origin coordinates, second clicks saves
        end coordinates, draws the line and saves its coordinates in a list."""

        # Checks if matplotlib is not using any tool (zoom or pan)
        tb = plt.get_current_fig_manager().toolbar
        if tb.mode == '' and event.inaxes == self.ax:
            self.xs.append(event.xdata)
            self.ys.append(event.ydata)
            if not self.clicked_once:
                marker, = self.ax.plot(self.xs, self.ys, color="orange", marker='.')
                marker.figure.canvas.draw()
                self.clicked_once = True
            else:
                line, = self.ax.plot(self.xs, self.ys, color="orange", marker='.', linewidth=1)
                line.figure.canvas.draw()
                self.linelist.append(line.get_data())
                self.fig.canvas.draw_idle()
                print(line.get_data())
                self.xs, self.ys = ([], [])
                self.clicked_once = False

    def get_values(self, event):
        """Returns the list of line coordinates when matplotlib window is closed."""
        return self.linelist

    def remove_selection(self, event):
        """Removes previous selection when r is pressed on keyboard"""
        if event.key == 'r':
            self.linelist.pop()
            print('Last line removed')

            self.ax.lines.pop()
            self.ax.lines.pop()


def beadlines(data):
    """Displays BeadlinesFigure with image from chosen path and returns values from get_values on shut down."""
    fig = BeadlinesFigure(data)
    return fig.get_values(fig)


##When the code is run, a navigator window will appear. Browse to and select folder containing .msr images to convert
print("Press Alt+tab to find the windows browser \n")
filename = easygui.diropenbox(default=os.path.expanduser("~Desktop"))
print(filename)

path=os.path.join(filename,"*.tif")

images = glob.glob(path)
print('There are ',len(images), 'Images in this folder')




#with MSRReader() as msrreader:

mapfoldernames={}

outpath=os.path.join(filename,"graphs")
os.makedirs(outpath, exist_ok=True)


for imagei in images:
    print(imagei)
    imagecomp = tifffile.imread(imagei)
    print(imagecomp.shape)
    
   
            #imagecomp=numpy.dstack(image1)


    imagecomp = numpy.moveaxis(imagecomp, 0, -1)
    print(imagecomp.shape)

# 
    
    nchannels=imagecomp.shape[2]
    print('nchannels',nchannels)
    if len(dispchannel)==1:
        img=imagecomp[:,:,dispchannel[0]]
        print('img',img.shape)
        linelist = beadlines(img)
    else:
        img=imagecomp[:,:,dispchannel]
        print('img',img.shape)
        img= numpy.moveaxis(img, -1, 0)
        print('img',img.shape)
        img = make_composite(img, luts,ranges=rangelut)
        linelist = beadlines(img)

    num_lines = len(linelist)
    print('num_lines',num_lines)


    figim, axesim = pyplot.subplots()
    axesim.axis('off')
    axesim.imshow(img, cmap='hot')
    with open(os.path.join(outpath,os.path.basename(imagei) +'_linelist.csv'),'w') as data:
        writerlines = csv.writer(data)
        for line in linelist:
            writerlines.writerow(line)
    for i, line in enumerate(linelist):
        fig, axes = pyplot.subplots(nrows=int(nchannels/2)+1, ncols=2, sharex=True,figsize=(6,12))
        axesim.plot(line[0], line[1],label='Line{}'.format(i))

        xcouple = (line[1][1], line[1][0])
        ycouple = (line[0][1], line[0][0])
        firstx = min(xcouple)
        firsty = ycouple[xcouple.index(firstx)]
        secondx = max(xcouple)
        secondy = ycouple[xcouple.index(secondx)]
        if xcouple[0]==xcouple[1]:
            secondx=xcouple[1]
            secondy=ycouple[1]
        print('secondx',secondx)
        print('secondy',secondy)
        traces = []
        
        for n in range(nchannels):
            trace = measure.profile_line(imagecomp[:,:,n], (firstx, firsty), (secondx, secondy), linewidth=Width)
            print(trace.shape)
            traces.append(trace)
        traces=numpy.array(traces)
        print(traces.shape)
        traces=traces/numpy.max(traces)
        xvalues=numpy.linspace(0,trace.shape[0], num=trace.shape[0])*pixelsize
        print(xvalues.shape)
        for n,trace in enumerate(traces):
            axes[int(nchannels/2),columns[n]].plot(xvalues,(trace/numpy.max(trace)),c=plotcolors[n], label='Line{} {}'.format(i,labels[n]))
            axes[rows[n],columns[n]].plot(xvalues,(trace/numpy.max(trace)),c=plotcolors[n], label='Line{}'.format(i))
            axes[rows[n],columns[n]].set_title(labels[n])
            axes[rows[n],columns[n]].set_xlim([0, 500])
            axes[rows[n],columns[n]].set_ylim([0, 1])
        axes[int(nchannels/2),0].set_xlim([0, 500])
        axes[int(nchannels/2),0].set_ylim([0, 1])
        axes[int(nchannels/2),1].set_xlim([0, 500])
        axes[int(nchannels/2),1].set_ylim([0, 1])
        axesim.legend()
        figim.savefig(os.path.join(outpath,os.path.basename(imagei)+'Image_Lines.pdf'.format(i)),transparent=True, bbox_inches='tight')
        fig.savefig(os.path.join(outpath,os.path.basename(imagei)+'Line{}_IntensityProfile.pdf'.format(i)),transparent=True, bbox_inches='tight')
        pyplot.close(fig)
        pyplot.close(figim)
    


