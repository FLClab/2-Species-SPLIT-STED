
"""
 This script allows the user to draw lines on an image and extract the intensity profile of the channels in the image along the line.
   The user can draw multiple lines on the image and the script will generate a graph for each line.
 The script will also generate an image showing the position of the lines drawn on the image.
   The user can choose the channels to display for line tracing in the image.
"""


import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot
import numpy
from skimage import measure
import os
import glob
import csv
import tifffile
from tiffwrapper import make_composite
import easygui


# rows and columns in the figure to plot each channel's line profile and the color of the line
rows= (0,1,0,1,2,0,1,2) 
columns= (2,2,0,0,0,1,1,1)
plotcolors = ('black','black', 'gray','blue','gold','gray','blue','gold')

#  plot titles
labels=["Mix Confocal","Mix STED","2 species phasor-FLIM - Confocal","2 species phasor-FLIM - STED","2 species SPLIT-STED","2 species phasor-FLIM - Confocal","2 species phasor-FLIM - STED","2 species SPLIT-STED"]

#channels to be used to make the composite image to be displayed in interactive window.
luts=("Cyan Hot", "Magenta Hot")
rangelut=((0, 100), (0, 100))
dispchannel=[4,7]

# Width of the selection line in pixels and the pixel size in nm
Width=2
pixelsize=20

# Dictionary to normalize the intensity of the channels
normch={0:[0],1:[1],2:[2,5],3:[3,6],4:[4,7],5:[2,5],6:[3,6],7:[4,7]}


class BeadlinesFigure:
    """Class that displays an image in a matplotlib window and allows the user to draw line segments.
    Line coordinates are saved in a list.
    Takes input image data as numpy array as parameter."""

    def __init__(self, img_data):

        self.fig, self.ax = pyplot.subplots(figsize=(15, 8))
        self.clicked_once = False
        self.xs = []
        self.ys = []
        self.linelist = []
        self.current_marker = None
        self.img_data = img_data

        pyplot.imshow(img_data, cmap='hot')
        #pyplot.colorbar()

        self.ax.set_title('Draw lines by clicking on endpoints. r removes last line ')

        self.fig.canvas.mpl_connect('button_press_event', self.click_event)
        self.fig.canvas.mpl_connect('key_press_event', self.remove_selection)
        self.fig.canvas.mpl_connect('close_event', self.get_values)

        pyplot.show()

    def click_event(self, event):
        """Allows the drawing of lines on the figure. First click saves origin coordinates, second clicks saves
        end coordinates, draws the line and saves its coordinates in a list."""

        # Checks if matplotlib is not using any tool (zoom or pan)
        tb = pyplot.get_current_fig_manager().toolbar
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


##When the code is run, a navigator window will appear. Browse to and select folder containing the tif images
print("Press Alt+tab to find the windows browser \n")
filename = easygui.diropenbox(default=os.path.expanduser("~Desktop"))
print(filename)

path=os.path.join(filename,"*.tif")

images = glob.glob(path)
print('There are ',len(images), 'Images in this folder')

#Create a folder to save the graphs
mapfoldernames={}
outpath=os.path.join(filename,"graphs")
os.makedirs(outpath, exist_ok=True)


for imagei in images:
    print(imagei)
    imagecomp = tifffile.imread(imagei)
    imagecomp = numpy.moveaxis(imagecomp, 0, -1)
    nchannels=imagecomp.shape[2]
# Display the image and ask the user to draw lines on the image from which to extract the intensity profile.
#If more than one channel is listed in dispchannel, the composite image will be displayed in the interactive window.
    if len(dispchannel)==1:
        img=imagecomp[:,:,dispchannel[0]]
        print('img',img.shape)
        linelist = beadlines(img)
#
    else:
        img=imagecomp[:,:,dispchannel]
        print('img',img.shape)
        img= numpy.moveaxis(img, -1, 0)
        print('img',img.shape)
        img = make_composite(img, luts,ranges=rangelut)
        linelist = beadlines(img)

    num_lines = len(linelist)
    print('num_lines',num_lines)

    # Display the image with the lines drawn on it
    figim, axesim = pyplot.subplots()
    axesim.axis('off')
    axesim.imshow(img, cmap='hot')

    # Save the line coordinates in a csv file
    with open(os.path.join(outpath,os.path.basename(imagei) +'_linelist.csv'),'w') as data:
        writerlines = csv.writer(data)
        for line in linelist:
            writerlines.writerow(line)
    for i, line in enumerate(linelist):
        fig, axes = pyplot.subplots(nrows=int((nchannels-2)/2)+1, ncols=3, sharex=True,figsize=(8,12))
        axesim.plot(line[0], line[1],label='Line{}'.format(i))

        xcouple = (line[1][1], line[1][0])
        ycouple = (line[0][1], line[0][0])
        firstx = min(xcouple)
        firsty = ycouple[xcouple.index(firstx)]
        secondx = max(xcouple)
        secondy = ycouple[xcouple.index(secondx)]
        # If the line is vertical, make the second point the one with the highest y value
        if xcouple[0]==xcouple[1]:
            secondx=xcouple[1]
            secondy=ycouple[1]
        print('secondx',secondx)
        print('secondy',secondy)
        traces = []
        # Extract the intensity profile of the channels along the line
        for n in range(nchannels):
            trace = measure.profile_line(imagecomp[:,:,n], (firstx, firsty), (secondx, secondy), linewidth=Width)
            print(trace.shape)
            traces.append(trace)
        traces=numpy.array(traces)
        print(traces.shape)
        # Normalize the intensity of the channels according to the dictionary normch
        for chnorm in normch.keys():
            chlist=normch[chnorm]
            print(chlist)
            print(traces[chlist].shape)
            traces[chnorm]=traces[chnorm]/numpy.max(traces[chlist])
      
        xvalues=numpy.linspace(0,trace.shape[0], num=trace.shape[0])*pixelsize
        print(xvalues.shape)
        for n,trace in enumerate(traces):

            axes[int((nchannels-2)/2),columns[n]].plot(xvalues,trace,c=plotcolors[n], label='Line{} {}'.format(i,labels[n]))
            axes[rows[n],columns[n]].plot(xvalues,trace,c=plotcolors[n], label='Line{}'.format(i))
            axes[rows[n],columns[n]].set_title(labels[n])
            #axes[rows[n],columns[n]].set_xlim([0, 500])
            axes[rows[n],columns[n]].set_ylim([0, 1])
        #axes[int(nchannels/2),0].set_xlim([0, 500])
        axes[int((nchannels-2)/2),0].set_ylim([0, 1])
        #axes[int(nchannels/2),1].set_xlim([0, 500])
        axes[int((nchannels-2)/2),1].set_ylim([0, 1])
        axesim.legend()
        figim.savefig(os.path.join(outpath,os.path.basename(imagei)+'Image_Lines.pdf'.format(i)),transparent=True, bbox_inches='tight')
        fig.savefig(os.path.join(outpath,os.path.basename(imagei)+'Line{}_IntensityProfile.pdf'.format(i)),transparent=True, bbox_inches='tight')
        pyplot.close(fig)
        pyplot.close(figim)
    


