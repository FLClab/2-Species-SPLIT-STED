import numpy
import os
import tifffile

try:
    from specpy import File
except:
    print("specpy not found, will use tifffile instead")
    pass


def select_channel(file,channel):
    if type(file)==dict:
        try:
            return file[channel]
        except KeyError:
            print(f"Channel {channel} not found in file.\n Available channels: {list(file.keys())}")
            return None
    else:
        #Selects a channel from a multi-channel image file.
        if len(file.shape) == 4:
            file=file[:,channel,:,:]  
            print(file.shape)
            return numpy.moveaxis(file, 0, -1)
        
        elif len(file.shape) == 3:
            print(file.shape)
            return numpy.moveaxis(file, 0, -1)
        
def load_image(file):
    if ".msr" in file:
        image=load_msr(file)
    else:
        image = tifffile.imread(file)
    return image

def load_msr(msrPATH):
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

def to_polar_coord(g, s, phi_IRF=0.4695955819269703, m_IRF=0.9527011687260826) :
    """
    Converts Cartesian phasor coordinates in the first harmonic to polar coordinates and performs calibration based on IRF measurement

    Parameters:
    g (list): List of x-coordinates.
    s (list): List of y-coordinates.
    phi_IRF (float): Mean angle of the IRF measurement in radians. Default is 0.4695955819269703.
    m_IRF (float): Mean magnitude of the IRF measurement. Default is 0.9527011687260826.

    Returns:
    tuple: A tuple containing two lists - m (magnitude) and phi (angle in radians).
    """

    m, phi = [], []
    for g_ij, s_ij in zip(g, s) :
        m.append(numpy.sqrt(g_ij ** 2 + s_ij ** 2))
        phi.append(numpy.arctan2(s_ij , g_ij))

    phi = [numpy.abs(i - phi_IRF) for i in phi] 
    m = [i / m_IRF for i in m]

    return m, phi

def to_polar_coord_2ndharm(g, s, phi_IRF=0.8874357044014783, m_IRF=0.9460309688519306) :
    """
    Converts Cartesian phasor coordinates in the first harmonic to polar coordinates and performs calibration based on IRF measurement

    Parameters:
    g (list): List of x-coordinates.
    s (list): List of y-coordinates.
    phi_IRF (float): Mean angle of the IRF measurement in radians. Default is 0.8874357044014783.
    m_IRF (float): Mean magnitude of the IRF measurement. Default is 0.9460309688519306.

    Returns:
    tuple: A tuple containing two lists - m (magnitude) and phi (angle in radians).
    """
    m, phi = [], []
    for g_ij, s_ij in zip(g, s) :
        m.append(numpy.sqrt(g_ij ** 2 + s_ij ** 2))
        phi.append(numpy.arctan2(s_ij , g_ij))
    phi = [i - phi_IRF for i in phi]
    m = [i / m_IRF for i in m] 

    return m, phi


def polar_to_cart(m, phi) :
    """
    Convert polar coordinates to Cartesian coordinates.

    Parameters:
    m (list): Magnitude values.
    phi (list): Angle values in radians.

    Returns:
    tuple: A tuple containing two lists - g (x-coordinate values) and s (y-coordinate values).
    """

    g, s = [], []
    for m_ij, phi_ij in zip(m, phi) :
        g.append(m_ij * numpy.cos(phi_ij))
        s.append(m_ij * numpy.sin(phi_ij))
    return g, s



from skimage import (
    filters, measure, morphology, segmentation)
from skimage.morphology import disk
#from skimage.morphology import (erosion, dilation, opening, closing,  # noqa
                               # white_tophat)

def get_foreground(image1):
    """
    Calculates the foreground of an image using triangle thresholding.

    Parameters:
    image1 (numpy.ndarray): The input image.

    Returns:
    numpy.ndarray: The thresholded image.

    """

    if len(image1.shape) == 3 :
        imsum = image1.sum(axis=2)
    else :
        imsum=image1
    #sigma =1
    #blurred = filters.gaussian(imsum, sigma=sigma)
    #blurred /= blurred.max()
    #blurred *= imsum
    seuil= filters.threshold_triangle(imsum)
    thresh =imsum > seuil
    #skimage.morphology.remove_small_objects(ar, min_size=64, connectivity=1, in_place=False, *, out=None)
    #footprint = morphology.disk(3)
    #dilated = morphology.dilation(thresh, footprint) * imsum
    return seuil



from numpy import ones,vstack
from numpy.linalg import lstsq

def line_equation(P1, P2):
    """
    Calculates the equation of a line given two points.

    Parameters:
    P1 (array-like): The coordinates of the first point (x1, y1).
    P2 (array-like): The coordinates of the second point (x2, y2).

    Returns:
    tuple: A tuple containing the slope (m) and y-intercept (c) of the line.

    Example:
    >>> P1 = [1, 2]
    >>> P2 = [3, 4]
    >>> line_equation(P1, P2)
    (1.0, 1.0)
    """
    A = vstack([P1, ones(len(P1))]).T
    m, c = lstsq(A, P2,rcond=None)[0]
    #print("y = {m}x + {c}".format(m=numpy.round(m,2),c=numpy.round(c,2)))
    return m, c

def ask_user(question):
    """Ask user a question and return the answer.

    Args:
        question (string): The question to ask the user.

    Returns:
        Bool: True if the answer is yes, False otherwise.
    """
    answer = input(question + " (Y/n) : ")
    if answer == "Y" or answer == "y" or answer == "yes" or answer == "Yes" or answer == "YES" or answer == "":
        return True
    else:
        return False