
"""
Script that performs 2 species SPLIT-STED phasor analysis on a set of images.
The script takes as input 2 folders of control images, one for each species, and a folder of mixed images.
It calculates the phasor distribution for each control image, builds a triangle in phasor space
 and then uses it to unmix the mixed images. 
 
 The script outputs:
    - the phasor distribution of the mixed images colored by the unmixing results
    - the mixed images colored by the unmixing results and the separated fractions
    - the phasor distribution of the control images
    - Composite image of the unmixed fluorophores with the third fraction removed
"""

import os
import glob
import numpy
import easygui
import pandas as pd
import time
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches
import os
import easygui
import seaborn
from sklearn.cluster import KMeans
import tifffile
import argparse
from sys import path as path1;
Functionspath=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Functions")
path1.append(Functionspath)
from Main_functions import (line_equation, to_polar_coord, polar_to_cart, load_image,select_channel, get_foreground)
from Phasor_functions import Median_Phasor,unmix3species
from tiffwrapper import make_composite,imsave,LifetimeOverlayer
from skspatial.objects import Circle
from skspatial.objects import Line
from skspatial.plotting import plot_2d


matplotlib.rcParams['axes.linewidth'] = 0.8


# ------------------ Default Input variables ----------------
params_dict = {
    "smooth_factor": 0.2,  # % 'Smoothing factor:'
    "phasor_smooth_cycles": 1,  # % 'Smoothing cycles phasor:'
    "foreground_threshold": 10,
    "harm1": 1,
}
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    name="nice-prism",
    colors=["#5F4690", "#1D6996", "#38A6A5", "#0F8554", "#73AF48", "#EDAD08", "#E17C05", "#CC503E"]
)
matplotlib.colormaps.register(cmap=cmap, force=True)
matplotlib.colormaps.register(cmap=cmap.reversed(), force=True)

# Define a custom argument type for a list of integers
def list_of_ints(arg):
    return list(map(int, arg.split(',')))

parser = argparse.ArgumentParser(description='Perform 2 species SPLIT-STED phasor analysis on a set of images')
parser.add_argument("-f1","--f1",help='Path to first folder', type=str )

parser.add_argument("-f2","--f2",help='Path to second folder',  type=str )
parser.add_argument("-f3","--f3",help='Path to second folder',  type=str )
parser.add_argument("-savefolder","--savefolder",help='Name of folder to save to',  type=str )
parser.add_argument("-numim","--numim",help='ID of images to use for controls',  type=list_of_ints )
args =vars(parser.parse_args())
print(args)
f1=args["f1"]
f2=args["f2"]
f3=args["f3"]
savefolder=args["savefolder"]
numimlist=args["numim"]

if None in [f1,f2,f3,savefolder]:

    # -----------------------------------------------------------


    # Select the folders containing the control images and the mixed images
    f1=easygui.diropenbox(default=os.path.expanduser("~Desktop"), title="Select folder containing the control images for the first fluorophore")
    f2=easygui.diropenbox(default=os.path.expanduser("~Desktop"), title="Select folder containing the control images for the second fluorophore")
    f3=easygui.diropenbox(default=os.path.expanduser("~Desktop"), title="Select folder containing the mixed images (the mixture of the two fluorophores)")

    savefolder=str(input("Name of Output folder: "))

    # Select the images to use for the controls. If none are selected, the script will ask for the image number

    #numimlist=[3,3,5,2,16,16,1,4]
    #numimlist=[27,27,12,3,6,6,9,3]
    #numimlist=[34,34,37,51,51,51,48,45]
    #numimlist=[6,6,1,4,5,5,3,1]
    #numimlist=[1,1,5,4,17,17,7,0] # PSD-Bassoon Cy5
    #numimlist=[0,0,1,9,22,22,7,0] # Spectrin Bassoon Cy5
    #numimlist=[15,15,0,5,22,22,7,0] #Tubulin Bassoon Cy5
    #numimlist=[7,7,1,0,0,0,1,19] # PSD Bassoon Cy3
    #numimlist = [1,1,0,8,8,8,9,7] #Actin Bassoon CY3
    #numimlist=[18,18,19,2,19,19,0,18] # Spectrin Bassoon Cy3
    #numimlist=[0,0,6,2,5,5,6,0] # Homer Bassoon Cy3 MediumAcq
    #numimlist=[0,0,5,1,4,4,3,2] # Homer Bassoon Cy3 ShortAcq
    #numimlist=[2,2,0,3,7,7,6,0] # Homer Bassoon Cy3 LongAcq
    

colors=['lightskyblue', 'deepskyblue','blue','mediumblue',"pink","lightpink",'lightcoral','indianred','springgreen']

labels=["Confocal","STED 10%","STED 20%","STED 30%","Confocal","STED 10%","STED 20%","STED 30%",'Mixture STED']
#labels=["Confocal","STED 5%","STED 10%","STED 15%","Confocal","STED 5%","STED 10%","STED 15%",'Mixture STED']

filenamescontrol = [f1,f1,f1,f1, f2,f2,f2,f2]
filenamemixed=f3
keys = [ 'Conf_635P {2}', 'STED_635P {2}', 'STED_635P {2}', 'STED_635P {2}',  'Conf_635P {2}', 'STED_635P {2}', 'STED_635P {2}', 'STED_635P {2}', 'STED_635P {2}']
#keys = ['Conf_ 594 {2}', 'STED_594 {2}', 'STED_594 {2}', 'STED_594 {2}','Conf_ 594 {2}', 'STED_594 {2}', 'STED_594 {2}', 'STED_594 {2}','STED_594 {2}']
#keys = ['Confocal_561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'Confocal_561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'STED 561 {11}']
keys=[0,1,1,1,0,1,1,1,1]
#keys = ['Conf_ 594 {2}', 'STED_594 {2}','STED_594 {2}','STED_594 {2}','Conf_ 594 {2}', 'STED_594 {2}','STED_594 {2}','STED_594 {2}']
#keys = ['Confocal_561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'Confocal_561 {11}', 'STED 561 {11}', 'STED 561 {11}', 'STED 561 {11}',  'STED 561 {11}']
#keys = [ 'Conf640 {10}', 'STED640 {10}', 'STED640 {10}', 'STED640 {10}','Conf640 {10}', 'STED640 {10}', 'STED640 {10}', 'STED640 {10}', 'STED640 {10}']
msrfiles = []
#plt.style.use('dark_background')
# -----------------------------------------------------------
#     Choisir un fichier msr parmi plusieurs

savefolder = os.path.join(os.path.expanduser("~/Desktop"), "Unmixing_"+savefolder+"_3Species")
os.makedirs(savefolder, exist_ok=True)


def abc_to_rgb(A=0.0,B=0.0,C=0.0):
# Map values A, B, C (all in domain [0,1]) to
# suitable red, green, blue values.
    return (min(B+C,1.0),min(A+C,1.0),min(A+B,1.0))
def plot_legend():
    ''' Plots a legend for the colour scheme
    given by abc_to_rgb. Includes some code adapted
    from http://stackoverflow.com/a/6076050/637562'''

    # Basis vectors for triangle
    basis = numpy.array([[0.0, 1.0], [-1.5/numpy.sqrt(3), -0.5],[1.5/numpy.sqrt(3), -0.5]])

    fig = plt.figure()
    ax = fig.add_subplot(111,aspect='equal')

    # Plot points
    a, b, c = numpy.mgrid[0.0:1.0:50j, 0.0:1.0:50j, 0.0:1.0:50j]
    a, b, c = a.flatten(), b.flatten(), c.flatten()

    abc = numpy.dstack((a,b,c))[0]
    #abc = filter(lambda x: x[0]+x[1]+x[2]==1, abc) # remove points outside triangle
    abc = list(map(lambda x: x/sum(x), abc)) # or just make sure points lie inside triangle ...

    data = numpy.dot(abc, basis)
    colours = [abc_to_rgb(A=point[0],B=point[1],C=point[2]) for point in abc]

    ax.scatter(data[:,0], data[:,1],marker=',',edgecolors='none',facecolors=colours,rasterized=True)

    # Plot triangle
    ax.plot([basis[_,0] for _ in range(3)],[basis[_,1] for _ in range(3)],**{'color':'black','linewidth':3})

    # Plot labels at vertices
    offset = 0.25
    fontsize = 32
    ax.text(basis[0,0]*(1+offset), basis[0,1]*(1+offset), '$A$', horizontalalignment='center',
            verticalalignment='center', fontsize=fontsize)
    ax.text(basis[1,0]*(1+offset), basis[1,1]*(1+offset), '$B$', horizontalalignment='center',
            verticalalignment='center', fontsize=fontsize)
    ax.text(basis[2,0]*(1+offset), basis[2,1]*(1+offset), '$C$', horizontalalignment='center',
            verticalalignment='center', fontsize=fontsize)    

    ax.set_frame_on(False)
    ax.set_xticks(())
    ax.set_yticks(())
    fig.savefig(os.path.join(savefolder,'Triangle_Legend.pdf'), transparent='True', bbox_inches="tight", dpi=900)
    fig.savefig(os.path.join(savefolder,'Triangle_Legend.png'), transparent='True', bbox_inches="tight", dpi=900)
    fig.savefig(os.path.join(savefolder,'Triangle_Legend.svg'), transparent='True', bbox_inches="tight", dpi=900)

plot_legend()

# Make list of all the images in the control folders and ask user to select the image to use as control
for k,filename in enumerate(filenamescontrol) :
    print(labels[k])
    
    extension = ".msr"
    path=os.path.join(filename,"*.msr")
    images = glob.glob(path)
    print('There are ',len(images), ' msr files in this folder')
    if len(images) == 0:
        path=os.path.join(filename,"*.tiff")
        images = glob.glob(path)
        print('There are ',len(images), ' tiff files in this folder')
        extension = ".tiff"

    for i,imagei in enumerate(images):
        print(i,os.path.basename(imagei)) 
    if numimlist is None:
        numim = int(input('Please enter the image number (1st=0): '))
        image = images[numim]
    else:
        image = images[numimlist[k]]
    msrfiles.append(image)
print(msrfiles)

path = os.path.join(filenamemixed, '*'+extension)
#path = os.path.join(filenamemixed, '*rep0*.msr')
mixedimages = glob.glob(path)

# Create figure for phasor plot
fig4,ax_scatter = plt.subplots(figsize=(2,2))
ax_scatter.set_xlim(0, 1.05)
ax_scatter.set_ylim(-0.05, 1)
ax_scatter.set_xlabel('g')
ax_scatter.set_ylabel('s')

# Draw the universal semi-circle
edge = numpy.linspace(start=0, stop=15, num=200)
theta = numpy.linspace(0, numpy.pi, 100)
r = 0.5
x1 = r*numpy.cos(theta) + 0.5
x2 = r*numpy.sin(theta)
ax_scatter.plot(x1,x2, color = "black", ls = "--",linewidth=0.8)


# Create a legend file
with open(os.path.join(savefolder,'legend.txt'),'w') as data:
    data.write("Controls\n")

# Loop through the control images and calculate their phasor distributions to build the reference triangle
scatterlist = []
barsxlist = []
barsylist = []
CoM_x, CoM_y = [], []
for i, msr in enumerate(msrfiles) : 
    df = pd.DataFrame(columns=['x','y'])
    dg = pd.DataFrame(columns=['g', 's'])
    imagemsr=load_image(msr)
    image1=select_channel(imagemsr,keys[i])
    print(image1.shape)
    # Write the control image name and the corresponding label and key to the legend file
    with open(os.path.join(savefolder,'legend.txt'),'a') as data:
        data.write("{}\t{}\t{}\n".format(labels[i],keys[i],msr))

    print(image1.shape)
    imsum = image1[:,:,10:111].sum(axis=2)
    imsum = imsum.astype('int16')
    
    seuil = get_foreground(imsum)


    print("Caclulation for an image of shape", image1.shape, "...")
    params_dict["foreground_threshold"] = seuil
    print("foreground_threshold=", params_dict["foreground_threshold"])
    # Calculate the phasor distribution of the foreground of the control image
    x,y,g_smoothed,s_smoothed, original_idxes= Median_Phasor(image1, params_dict, **params_dict)
    df['x']=x.flatten()
    df['y']=y.flatten()
    m, phi = to_polar_coord(df['x'], df['y'])
    g,s =polar_to_cart(m, phi)
    dg['g'], dg['s'] = g, s
    # Calculate the centroid of the phasor distribution using KMeans clustering
    kmeans = KMeans(n_clusters = 1, init = 'k-means++', random_state = 42)
    y_kmeans = kmeans.fit_predict(dg)
    # Store the centroid coordinates in a list
    CoM_x.extend(kmeans.cluster_centers_[:, 0][:].tolist())
    CoM_y.extend(kmeans.cluster_centers_[:, 1][:].tolist())
    # Plot the phasor distribution of the control image
    a=ax_scatter.scatter(g, s, s=0.5, c=colors[i], alpha=0.10,label=labels[i],rasterized=True)
    scatterlist.append(a)

# Find nearest points to the semi-circle of the confocal phasor centroids
xaxis = numpy.linspace(0, 1.5, 100)
norm = numpy.sqrt((CoM_x[0] - 0.5) ** 2 + (CoM_y[0] ** 2))
Pn_x = 0.5 + (r * (CoM_x[0] - 0.5) / norm)
Pn_y = 0 + r * (CoM_y[0] - 0) / norm
P_n = numpy.array([Pn_x, Pn_y])
norm = numpy.sqrt((CoM_x[4] - 0.5) ** 2 + (CoM_y[4] ** 2))
P2_x = 0.5 + (r * (CoM_x[4] - 0.5) / norm)
P2_y = 0 + r * (CoM_y[4] - 0) / norm
p2 = numpy.array([P2_x, P2_y])

# Separate the centroids of the two species into 2 lists
PointsSpecies1 = numpy.stack(
    [numpy.array([Pn_x, Pn_y]), numpy.array([CoM_x[1], CoM_y[1]]), numpy.array([CoM_x[2], CoM_y[2]]),
        numpy.array([CoM_x[3], CoM_y[3]])])
PointsSpecies2 = numpy.stack(
    [numpy.array([P2_x, P2_y]), numpy.array([CoM_x[5], CoM_y[5]]), numpy.array([CoM_x[6], CoM_y[6]]),
        numpy.array([CoM_x[7], CoM_y[7]])])


# Fit lines to the phasor centroids of each of the two species
coeffs1 = numpy.polyfit([Pn_x, CoM_x[1], CoM_x[2], CoM_x[3]], [Pn_y, CoM_y[1], CoM_y[2], CoM_y[3]], 1)
coeffs2 = numpy.polyfit([P2_x, CoM_x[5], CoM_x[6], CoM_x[7]], [P2_y, CoM_y[5], CoM_y[6], CoM_y[7]], 1)
y1 = coeffs1[0] * xaxis + coeffs1[1]
y2 = coeffs2[0] * xaxis + coeffs2[1]

# Find the intersection point of the two lines
det = coeffs2[0] - coeffs1[0]
x = (coeffs1[1] - coeffs2[1]) / det
y = (coeffs2[0] * coeffs1[1] - coeffs1[0] * coeffs2[1]) / det
p0 = numpy.array([x, y])
print('p0', p0)

# Check if the intersection point is inside the universal semi-circle
circ = Circle([0.5, 0], 0.5)
check = circ.contains_point([x, y])
if check == False:  # If intersection point is outside the circle, find intersection with circle
    print("I'm outside the circle, coming in!")
    circle = Circle([0.5, 0], 0.5)
    line1 = Line.from_points([Pn_x, Pn_y], [x, y])
    line2 = Line.from_points([P2_x, P2_y], [x, y])
    point_a, point_b = circle.intersect_line(line1)
    point_c, point_d = circle.intersect_line(line2)
    print('abcd', point_a, point_b, point_c, point_d)
    print(numpy.array([point_a, point_c]))
    p0 = numpy.mean(numpy.array([point_a, point_c]), axis=0)
    x = p0[0]
    y = p0[1]
if y < 0:  # If intersection point is under the semi-circle, find intersection with x-axis
    print("I'm in the negatives, coming up!")
    line1 = Line.from_points([Pn_x, Pn_y], [x, y])
    line2 = Line.from_points([P2_x, P2_y], [x, y])
    line3 = Line.from_points([0, 0], [1,0])
    point_a = line3.intersect_line(line1)
    point_b = line3.intersect_line(line2)
    p0 = numpy.mean(numpy.array([point_a, point_b]), axis=0)
    x = p0[0]
    y = p0[1]
print('p0', p0)
print("POINTS", P_n, p2, p0)

# Save the phasor containing the phasor distributions of the control images
fig4.savefig(os.path.join(savefolder, 'Phasor_SeparateSTED_3species_LineControls_controlsonly.pdf'),
                transparent='True', bbox_inches="tight",dpi=900)

# Clear the scatter plot
lines = [mpatches.Patch(color=colors[j], label=labels[j]) for j in range(len(labels))]
t = [scatter.remove() for scatter in scatterlist]

# Plot control image centroids, the reference triangle with the control points and lines and save the figure
pnscatter = ax_scatter.scatter(Pn_x, Pn_y, s=50, c='gold')
p2scatter = ax_scatter.scatter(P2_x, P2_y, s=50, c='gold')
p0scatter = ax_scatter.scatter(p0[0], p0[1], s=50, c='gold')
p0pnline = ax_scatter.plot([Pn_x, p0[0]], [Pn_y, p0[1]], c='dodgerblue')
p2pnline = ax_scatter.plot([Pn_x, P2_x], [Pn_y, P2_y], c='dodgerblue')
p0p2line = ax_scatter.plot([P2_x, p0[0]], [P2_y, p0[1]], c='dodgerblue')
centroidscatter1=ax_scatter.scatter(PointsSpecies1[:,0],PointsSpecies1[:,1],s=50,c=colors[0:4])
centroidscatter2 = ax_scatter.scatter(PointsSpecies2[:,0],PointsSpecies2[:,1], s=50,c=colors[4:-1])
fig4.savefig(os.path.join(savefolder, "Phasor_3species_LineControls_ControlsOnly.pdf"), transparent='True',
                bbox_inches="tight")
# Clear the scatter plot
pnscatter.remove()
p2scatter.remove()
p0scatter.remove()
centroidscatter1.remove()
centroidscatter2.remove()
ax_scatter.lines[-1].remove()
ax_scatter.lines[-1].remove()
ax_scatter.lines[-1].remove()


# Loop through the mixed images and unmix them using the phasor triangle built from the control images
for m,mixedimage in enumerate(mixedimages):

    d_melange = pd.DataFrame(columns=['g', 's'])
    df = pd.DataFrame(columns=['x', 'y'])
    dg = pd.DataFrame(columns=['g', 's'])
    print("***********************************************************")
    print("Working on image number ",m," out of ",len(mixedimages))
    print("***********************************************************")
    imagemsr = load_image(mixedimage)
    print(mixedimage)


    image1= select_channel(imagemsr,keys[-1])
    print(image1.shape)
    imsum = image1[:,:,10:111].sum(axis=2)
    imsum = imsum.astype('int16')

    seuil = 5
    params_dict["foreground_threshold"] = seuil
    print("foreground_threshold=", params_dict["foreground_threshold"])
    # Calculate the phasor distribution of the foreground of the mixed image
    x, y, g_smoothed, s_smoothed, original_idxes = Median_Phasor(image1, params_dict, **params_dict)
    df['x'] = x.flatten()
    df['y'] = y.flatten()
    m, phi = to_polar_coord(df['x'], df['y'])
    g, s = polar_to_cart(m, phi)
    dg['g'], dg['s'] = g, s
    d_melange['g'], d_melange['s'] = g, s
    # Calculate 2 centroids of the phasor distribution using KMeans clustering
    kmeans = KMeans(n_clusters=2, init='k-means++', random_state=42)
    y_kmeans = kmeans.fit_predict(dg)
    p3 = d_melange[['g', 's']].to_numpy() 
    print("p3",p3.shape)

    # Unmix the mixed image using the phasor triangle built from the control images
    imsum_flat_lin1, imsum_flat_lin2, imsum_flat_lin3, Solve = unmix3species(p3, original_idxes, image1, P_n, p2,p0)


    print('F1',imsum_flat_lin1.min(),imsum_flat_lin1.max())
    print('F2',imsum_flat_lin2.min(),imsum_flat_lin2.max())
    print('F3',imsum_flat_lin3.min(),imsum_flat_lin3.max())

    # Plot the phasor distribution of the mixed image colored by the unmixing results
    colours = [abc_to_rgb(A=point[0],B=point[1],C=point[2]) for point in numpy.transpose(Solve)]
    print("colours",numpy.min(colours),numpy.max(colours))
    mixphasor = ax_scatter.scatter(p3[:,0],p3[:,1],s=1,c=colours,rasterized=True)
    # Save the phasor distribution of the mixed image as a png
    fig4.savefig(os.path.join(savefolder, "Phasor_3species_LinesControls_{}.png".format(os.path.basename(mixedimage).split(extension)[0])), transparent='True',
                    bbox_inches="tight",dpi=900)

    # Plot the control points and the lines connecting them and save the figure as a pdf
    pnscatter = ax_scatter.scatter(Pn_x, Pn_y, s=50, c='gold')
    p2scatter = ax_scatter.scatter(P2_x, P2_y, s=50, c='gold')
    p0scatter = ax_scatter.scatter(p0[0], p0[1], s=50, c='gold')
    p0pnline = ax_scatter.plot([Pn_x, p0[0]], [Pn_y, p0[1]], c='dodgerblue')
    p2pnline = ax_scatter.plot([Pn_x, P2_x], [Pn_y, P2_y], c='dodgerblue')
    p0p2line = ax_scatter.plot([P2_x, p0[0]], [P2_y, p0[1]], c='dodgerblue')
    
    fig4.savefig(os.path.join(savefolder, "Phasor_3species_LinesControls_{}.pdf".format(os.path.basename(mixedimage).split(extension)[0])), transparent='True',
                    bbox_inches="tight",dpi=900)
    
    # Clear the scatter plot
    mixphasor.remove()
    pnscatter.remove()
    p2scatter.remove()
    p0scatter.remove()
    ax_scatter.lines[-1].remove()
    ax_scatter.lines[-1].remove()
    ax_scatter.lines[-1].remove()

    
    imsum_flat_lin3[numpy.isinf(imsum_flat_lin3)] = 0
    imsum_flat_lin3[numpy.isnan(imsum_flat_lin3)] = 0
    imsum_flat_lin3 = numpy.where(imsum_flat_lin3 < 0, 0, imsum_flat_lin3)
    imsum_flat_lin3 = numpy.where(imsum_flat_lin3 > 1, 1, imsum_flat_lin3)
    q = 1 - imsum_flat_lin3
    imsum_flat_bi = imsum_flat_lin2 / q
    imsum_flat_bi[numpy.isinf(imsum_flat_bi)] = 0
    imsum_flat_bi[numpy.isnan(imsum_flat_bi)] = 0
    imsum_flat_bi = numpy.where(imsum_flat_bi < 0, 0, imsum_flat_bi)
    imsum_flat_bi = numpy.where(imsum_flat_bi > 1, 1, imsum_flat_bi)

    imsum_flat_bi1 = imsum_flat_lin1 / q
    imsum_flat_bi1[numpy.isinf(imsum_flat_bi1)] = 0
    imsum_flat_bi1[numpy.isnan(imsum_flat_bi1)] = 0
    imsum_flat_bi1 = numpy.where(imsum_flat_bi1 < 0, 0, imsum_flat_bi1)
    imsum_flat_bi1 = numpy.where(imsum_flat_bi1 > 1, 1, imsum_flat_bi1)

    
    # Apply the calculated fractions to the intensity image to obtain the fraction images
    fraction1 = imsum_flat_bi.copy()
    fraction2 = imsum_flat_bi1.copy()
    fraction3=imsum_flat_lin3.copy()
    fraction3 *= imsum
    difference = imsum - fraction3
    fraction1 *= difference
    fraction2 *= difference


    # Plot the unmixed images and the intensity image
    fig_im3, ax_im3 = plt.subplots(ncols=5,nrows=1,figsize=(20,6))
    ax_im3[0].axis('off')
    ax_im3[1].axis('off')
    ax_im3[2].axis('off')
    ax_im3[3].axis('off')
    ax_im3[4].axis('off')
    ax_im3[0].set_title('Intensity image', fontsize=16)
    ax_im3[1].set_title('Fraction 1', fontsize=16)
    ax_im3[2].set_title('Fraction 2', fontsize=16)
    ax_im3[3].set_title('Fraction 3', fontsize=16)
    ax_im3[4].set_title('Fraction 1 and 2 composite', fontsize=16)
    imsum_disp =ax_im3[0].imshow(imsum, cmap='hot')
    imsum_flat3 =ax_im3[1].imshow(fraction2, cmap='hot')
    imsum_flat5 =ax_im3[2].imshow(fraction1, cmap='hot')
    imsum_flat6 =ax_im3[3].imshow(fraction3, cmap='hot')

    cbar =fig_im3.colorbar(imsum_disp,ax=ax_im3[0], fraction=0.025, pad=0.01)
    #cbar.set_label("Intensity", fontsize=12)
    cbar =fig_im3.colorbar(imsum_flat3,ax=ax_im3[1], fraction=0.025, pad=0.01)
    #cbar.set_label("Intensity", fontsize=12)
    cbar2 =fig_im3.colorbar(imsum_flat5,ax=ax_im3[2], fraction=0.025, pad=0.01)
    #cbar2.set_label("Intensity", fontsize=12)
    cbar3 =fig_im3.colorbar(imsum_flat6,ax=ax_im3[3], fraction=0.025, pad=0.01)
    #cbar3.set_label("Intensity", fontsize=12)

    
    lifetime=numpy.dstack(( imsum_flat_lin1, imsum_flat_lin2, imsum_flat_lin3))


    red=lifetime[:,:,1]+lifetime[:,:,2]
    green=lifetime[:,:,0]+lifetime[:,:,2]
    blue=lifetime[:,:,0]+lifetime[:,:,1]
    
    lifetime_rgb =numpy.dstack((red,green,blue))

    print(lifetime_rgb.shape,lifetime_rgb.dtype)
    print(numpy.min(lifetime_rgb[:,:,0]),numpy.max(lifetime_rgb[:,:,0]))  
    print(numpy.min(lifetime_rgb[:,:,1]),numpy.max(lifetime_rgb[:,:,1]))  
    print(numpy.min(lifetime_rgb[:,:,2]),numpy.max(lifetime_rgb[:,:,2]))  

    plt.figure()
    plt.imshow(lifetime_rgb)
    plt.axis('off')
    plt.savefig(os.path.join(savefolder, os.path.basename(mixedimage).split(extension)[0] + "_STED3species_LineControls_lifetimergb.pdf"), transparent='True', bbox_inches="tight")

    lifetime_rgb=numpy.moveaxis(lifetime_rgb, 2, 0)
    print(lifetime_rgb.shape)
    filenameout = os.path.join(savefolder, os.path.basename(mixedimage).split(extension)[
        0] + "_STED3species_LineControls_lifetimergb.tiff")
    tifffile.imwrite(filenameout, lifetime_rgb.astype(numpy.float32), photometric='rgb')

    overlayer = LifetimeOverlayer(lifetime, imsum/imsum.max(), cname='coolwarm')
    lifetime_rgb = overlayer.get_overlay_RGB(
        lifetime_minmax=(0., 1),
        intensity_minmax=(0, 0.5) # inTensity saturated to get more bright regions
                )

    imsum_flattt =ax_im3[4].imshow(lifetime_rgb)

# Saved the unmixed fraction images as composite tiff images
    imagecomp = numpy.dstack((fraction1,fraction2))
    imagecomp = numpy.moveaxis(imagecomp, 2, 0)
    filenameout = os.path.join(savefolder, os.path.basename(mixedimage).split(extension)[
        0] + "_STED3species_LineControls_UnmixedComposite.tiff")
    imsave(file=filenameout, data=imagecomp.astype(numpy.uint16), composite=True, luts=("Magenta Hot", "Cyan Hot"),
            pixelsize=(20E-3, 20E-3))
    
# Save the intensity of the mixed image
    filenameout = os.path.join(savefolder,
                                os.path.basename(mixedimage).split(extension)[0] + "_STED3species_LineControls_MixedIntensity.tiff")
    print(filenameout)
    imsave(file=filenameout, data=imsum.astype(numpy.uint16), luts="Red Hot", pixelsize=(20E-3, 20E-3))

# Save the unmixed fractions as tiff images
    filenameout = os.path.join(savefolder,
                                os.path.basename(mixedimage).split(extension)[0] + "_STED3species_LineControls_f1f2f3.tiff")
    imagecomp = numpy.dstack((imsum_flat_bi, imsum_flat_bi1,imsum_flat_lin3))
    imagecomp = numpy.moveaxis(imagecomp, 2, 0)
    tifffile.imwrite(filenameout, imagecomp)


    filenameout = os.path.join(savefolder,
                                os.path.basename(mixedimage).split(extension)[0] + "_STED3species_LineControls_F1Overlay.tiff")

    tifffile.imwrite(filenameout, lifetime_rgb.astype(numpy.float32))

    fig_im3.savefig(os.path.join(savefolder, 'Images_SeparateSTED3species_LineControls_' +
                                    os.path.basename(mixedimage).split(extension)[0] + '.pdf'), transparent='True',
                    bbox_inches="tight")


