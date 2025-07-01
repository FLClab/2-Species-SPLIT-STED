""" Computes the Full-Width at half maximum (FWHM) of the Instrument Response Function (IRF)
measurement.  
    
    Use a Confocal-FLIM image of a gold bead sample that backscatters the laser to the detector.
    The script will plot the histogram of the IRF measurement and fit a Gaussian to it.
    The FWHM is calculated from the fitted Gaussian.

"""

import matplotlib.pyplot as plt
import glob
import numpy
import easygui
from scipy.optimize import curve_fit
from sklearn.cluster import KMeans
import pandas
import os.path
from sys import path as path1;
Functionspath=os.path.join(os.path.dirname(os.path.dirname(__file__)), "Functions")
path1.append(Functionspath)
from Main_functions import get_foreground,load_image,select_channel, to_polar_coord
from Phasor_functions import Median_Phasor
# ------------------ Default Input variables ----------------
params_dict = {
    # Parameter in option in the matlab code
    #    "Tg" : 6, #% 'First frame to sum:'
    "Nb_to_sum": 250,  # The Tg infered from this variable override Tg
    "smooth_factor": 2,  # % 'Smoothing factor:'
    "im_smooth_cycles": 0,  # % 'Smoothing cycles image:'
    "phasor_smooth_cycles": 1,  # % 'Smoothing cycles phasor:'
    "foreground_threshold": 10,
    "tau_exc": numpy.inf,  # % 'Tau_exc'
    "intercept_at_origin": False,  # % 'Fix Intercept at origin'

    # Parameters that are calculated in th matlab code but that could be modified manually
    "M0": None,
    "Min": None,

    # Paramaters that are fixed in the matlab code
    "m0": 1,
    "harm1": 1,  # MATLAB code: harm1=1+2*(h1-1), where h1=1
    "klog": 4,
}

# -----------------------------------------------------------
# Select the folder containing the images
filename=os.path.join(os.path.dirname(__file__), "IRF_measurement")
key= 'STAR 635P_CONF {0}'
#key=0

# Make list of all the images in the folder
extension = ".msr"
path=os.path.join(filename,"*.msr")
images = glob.glob(path)
print('There are ',len(images), ' msr files in this folder')
if len(images) == 0:
    path=os.path.join(filename,"*.tiff")
    images = glob.glob(path)
    print('There are ',len(images), ' tiff files in this folder')
    extension = ".tiff"



# -----------------------------------------------------------
#   Ask user to select the image file to extract the IRF

for imagei in images:
    print(os.path.basename(imagei))
# Ask user to choose an image file
numim = int(input('Enter the index of the file to load (1st=0): '))
images = images[numim]
# Read the selected image file and extract the correct channel image
imagemsr = load_image(images)
image1=select_channel(imagemsr, key)
dim = image1.shape
# Find foreground threshold using the Triangle method

imsum= numpy.sum(image1, axis=2)
seuil = get_foreground(imsum)

# Plot the image and the foreground mask
fig, ax = plt.subplots(ncols=2,nrows=1, figsize=(10, 5))
imgplot1 = ax[0].imshow(imsum, cmap='hot')
cbar =fig.colorbar(imgplot1)
ax[1].imshow(imsum > seuil) 
ax[0].axis('off')
ax[1].axis('off')


# -----------------------------------------------------------
# Calculate a gaussian fit to the histogram of the IRF measurement's foreground
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
# Plot the histogram and the fitted Gaussian
y /= y.max()
ax_gaus.plot(absci, y , 'b+:', label='DATA')
ax_gaus.plot(absci,gaussian,'ro:',label='fit')
ax_gaus.set_xlabel('time [ns]', fontsize = 20)
ax_gaus.set_ylabel('Normalised intensity', fontsize = 20)
ax_gaus.xaxis.set_tick_params(labelsize=20)
ax_gaus.yaxis.set_tick_params(labelsize=20)

ax_gaus.legend(fontsize = 20)

#Create the figure and axes for the phasor plot
fig_centroids,ax_centroids = plt.subplots(figsize=(4,4))
theta = numpy.linspace(0, numpy.pi, 100)
r = 0.5
x1 = r * numpy.cos(theta) + 0.5
x2 = r * numpy.sin(theta)
ax_centroids.plot(x1, x2, color="black", ls="--",linewidth=0.8)
ax_centroids.set_xlim(-0.02, 1.02)
ax_centroids.set_ylim(-0.02, 1.02)
ax_centroids.set_xlabel('g')
ax_centroids.set_ylabel('s')


# Calculate the phasor distribution for all pixels in the image that are above the foreground threshold
df = pandas.DataFrame(columns=['x','y'])
dg = pandas.DataFrame(columns=['g', 's'])
print("Caclulation for an image of shape", image1.shape, "...")

params_dict["foreground_threshold"]=20
print("foreground_threshold=", params_dict["foreground_threshold"])
x,y,g_smoothed,s_smoothed, orginal_idxs= Median_Phasor(image1, params_dict, **params_dict, show_plots=False)
df['x']=x.flatten()
df['y']=y.flatten()
ax_centroids.scatter(df['x'],df['y'],s=0.5, c="gray",alpha=0.1,rasterized=True)

# Calculate the centroid of the phasor distribution
kmeans = KMeans(n_clusters=1, init='k-means++', random_state=42)
y_kmeans = kmeans.fit_predict(df)
CoM_x=kmeans.cluster_centers_[:, 0][:]
CoM_y=kmeans.cluster_centers_[:, 1][:]

ax_centroids.scatter(CoM_x, CoM_y, color="k", s=25)

# Convert the centroid coordinates to polar coordinates
m, phi = to_polar_coord(CoM_x, CoM_y,phi_IRF=0,m_IRF=1)

# Print the centroid coordinates
print("Centroid coordinates (g,s):", CoM_x, CoM_y)
print("Centroid coordinates in polar coordinates (m,phi):", m, phi)

plt.show()

