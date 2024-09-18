import numpy
import scipy

from sklearn.cluster import KMeans
import skimage.io as skio
import skimage
import plotly.express as px
from plotly.offline import plot
import dtcwt
import matplotlib.pyplot as plt
params_dict = {
    # Parameter in option in the matlab code
        "Tg" : 10, #% 'First frame to sum:'
    "Nb_to_sum": 250,  # The Tg infered from this variable override Tg
    "smooth_factor": 0.2,  # % 'Smoothing factor:'
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
def Median_Filter(im, sm, n):
    """
    Applies a median filter to the input image.

    Parameters:
    - im: numpy.ndarray
        The input image to be filtered.
    - sm: float
        The smoothing factor for the filter.
    - n: int
        The number of times the filter should be applied.

    Returns:
    - numpy.ndarray
        The filtered image.
    """
    #Add padding to get rid of border artifact. Padding is removed before returning result
    im=numpy.pad(im,3,mode='edge')
    filt = 1 / (8 + 1 / sm) * numpy.array([[1, 1, 1],
                                        [1, 1 / sm, 1],
                                        [1, 1, 1], ])
    for _ in range(n):
        im = scipy.signal.convolve2d(im, filt, mode='same')
    return im[3:-3,3:-3]


def Median_Phasor(sted_stack_fname, params_dict, Nb_to_sum, smooth_factor, foreground_threshold,
               im_smooth_cycles, phasor_smooth_cycles, tau_exc, intercept_at_origin, m0, M0, Min, harm1,
               klog, show_plots=False, new_tau_exc=None, return_analysis_fig=False):
    """
    Calculate the phasor distribution for a given FLIM image with median filtering.

    Args:
        sted_stack_fname (str or ndarray): The filename or the STED stack itself.
        params_dict (dict): A dictionary containing various parameters.
        Nb_to_sum (int): The number of frames to sum.
        smooth_factor (int): The smoothing factor.
        foreground_threshold (float): The threshold for foreground detection.
        im_smooth_cycles (int): The number of smoothing cycles for the image.
        phasor_smooth_cycles (int): The number of smoothing cycles for the phasor.
        tau_exc (float): The excitation lifetime.
        intercept_at_origin (bool): Whether to intercept at the origin.
        m0 (float): The m0 value.
        M0 (float): The M0 value.
        Min (float): The Min value.
        harm1 (int): The harmonic value.
        klog (float): The klog value.
        show_plots (bool, optional): Whether to show plots. Defaults to False.
        new_tau_exc (float, optional): The new excitation lifetime. Defaults to None.
        return_analysis_fig (bool, optional): Whether to return the analysis figure. Defaults to False.

    Returns:
        tuple: A tuple containing x, y, g_smoothed, s_smoothed, and original_idxes.
    """
    
    if new_tau_exc:
        tau_exc = new_tau_exc

    if type(sted_stack_fname) == str:
        sted_stack = skio.imread(sted_stack_fname)
        sted_stack = numpy.moveaxis(sted_stack, 0, -1)
    else:
        sted_stack = sted_stack_fname

    X, Y, N = sted_stack.shape
    if Nb_to_sum:
        Tg = N - Nb_to_sum + 1

    g = (sted_stack * numpy.cos(2 * numpy.pi * (harm1) * numpy.arange(N) / N)).sum(axis=2) / sted_stack.sum(axis=2)
    g[numpy.isnan(g)] = 0
    s = (sted_stack * numpy.sin(2 * numpy.pi * (harm1) * numpy.arange(N) / N)).sum(axis=2) / sted_stack.sum(axis=2)
    s[numpy.isnan(s)] = 0

    g = Median_Filter(g, smooth_factor, phasor_smooth_cycles)
    s = Median_Filter(s, smooth_factor, phasor_smooth_cycles)

    modules_smoothed = numpy.sqrt(s ** 2 + g ** 2)

    Ngat = sted_stack[:, :, 10:111].sum(axis=2)
  
    modules_smoothed[Ngat < foreground_threshold] = 0

    original_idxes = numpy.arange(len(g.flatten()))
    original_idxes = original_idxes[numpy.logical_and(Ngat > foreground_threshold, modules_smoothed > 0).flatten()]

    x = g[numpy.logical_and(Ngat > foreground_threshold, modules_smoothed > 0)]
    y = s[numpy.logical_and(Ngat > foreground_threshold, modules_smoothed > 0)]
    g_smoothed = g[:, :]
    s_smoothed = s[:, :]
    return x,y,g_smoothed,s_smoothed, original_idxes


def DTCWT_Phasor(sted_stack_fname, foreground_threshold,nlevels=2, neighborhood=50):
    """
    Compute the Phasor with Complex Wavelet filtering for a given FLIM image stack.
    Reproduced from https://doi.org/10.1364/BOE.420953

    Parameters:
    sted_stack_fname (str or ndarray): The filename of the FLIM image or the image itself.
    foreground_threshold (float): The threshold for foreground detection.
    nlevels (int): The number of levels for the DTCWT transform (default: 2).
    neighborhood (int): The size of the neighborhood for local noise estimation (default: 50).

    Returns:
    tuple: A tuple containing the G, S, original_idxes, Images, and Images_Filtered.

    - G (ndarray): The G component of the Phasor.
    - S (ndarray): The S component of the Phasor.
    - original_idxes (ndarray): The original indices of the Phasor in the image.
    - Images (list): A list of the original images used in the computation.
    - Images_Filtered (list): A list of the filtered images used in the computation.
    """
    
    if type(sted_stack_fname) == str:
        sted_stack = skio.imread(sted_stack_fname)
        sted_stack = numpy.moveaxis(sted_stack, 0, -1)
    else:
        sted_stack = sted_stack_fname
    sted_stack[numpy.sum(sted_stack[:,:,10:111],axis=2)<foreground_threshold,:]=0
    X, Y, N = sted_stack.shape
    Im=sted_stack.sum(axis=2)
    FReal = (sted_stack * numpy.cos((2 * numpy.pi * numpy.arange(N))/ N)).sum(axis=2)
    FImag = (sted_stack * numpy.sin((2 * numpy.pi * numpy.arange(N))/ N)).sum(axis=2)
    #print("Im",numpy.min(Im),numpy.max(Im))
    FReal[FReal <0]=0
    FImag[FImag < 0]=0

    Images=[Im,FReal,FImag]
    Images_Filtered=[]
    transform = dtcwt.Transform2d(biort='legall',qshift='qshift_a')
    for I in Images:

        IAscombe = 2 * numpy.sqrt(numpy.abs(I) + (3 / 8))
        IAscombe[numpy.isnan(IAscombe)] = 0
        I_t=transform.forward(IAscombe,nlevels=nlevels)

        lowpass=(I_t.lowpass)
        Levels=[I_t.highpasses[i] for i in range(nlevels)]
        Amplitudes=[numpy.absolute(Level) for Level in Levels]

        #GlobalNoise=numpy.median([Amplitudes[0][:,:,1],Amplitudes[0][:,:,4]])/0.6745
        GlobalNoise = numpy.median([Amplitudes[0][:, :, 0], Amplitudes[0][:, :, 2]]) / 0.6745


        LocalNoiseEstimate=[]
        for l in range(nlevels):
            localnoiseband=numpy.zeros(Amplitudes[l].shape)
            for b in range(6):
                currentAmp=Amplitudes[l][:,:,b]**2

                padded = numpy.pad(currentAmp, neighborhood, mode='edge')
                localnoise=1/(2*neighborhood+1)**2*numpy.sum(numpy.lib.stride_tricks.sliding_window_view(padded,(2*neighborhood+1,2*neighborhood+1)),axis=(-1,-2))

                #kernel = numpy.ones((2 * neighborhood + 1, 2 * neighborhood + 1))
                #localnoise=1/(2*neighborhood+1)**2*scipy.signal.convolve2d(kernel,padded,mode='valid')
                localnoiseband[:,:,b]=localnoise
            LocalNoiseEstimate.append(localnoiseband)
        FilteredPyramid=[]
        for l in range(nlevels-1):
            filteredcoeffsL=numpy.zeros(Amplitudes[l].shape,dtype=numpy.complex_)

            for b in range(6):
                currentLevel = Levels[l][:, :, b]
                currentAmp=Amplitudes[l][:,:,b]
                currentParentAmp=Amplitudes[l+1][:,:,b]
                currentParentAmpreps=numpy.repeat(currentParentAmp,2,axis=0)
                currentParentAmpreps = numpy.repeat(currentParentAmpreps, 2, axis=1)
                if currentAmp.shape[0]<currentParentAmpreps.shape[0]:
                    currentParentAmpreps=numpy.delete(currentParentAmpreps,-1,0)
                if currentAmp.shape[1]<currentParentAmpreps.shape[1]:
                    currentParentAmpreps=numpy.delete(currentParentAmpreps,-1,1)


                currentlocalnoise=LocalNoiseEstimate[l][:,:,b]
                denom=numpy.sqrt((currentAmp**2+currentParentAmpreps)*numpy.absolute(currentlocalnoise-GlobalNoise))
                filteredcoeffsxy=currentLevel*(1-((numpy.sqrt(3)*GlobalNoise**2)/denom))

                filteredcoeffsL[:,:,b]=filteredcoeffsxy
            FilteredPyramid.append(filteredcoeffsL)

        FilteredPyramid.append(Amplitudes[nlevels-1])
        pyramidd=dtcwt.Pyramid(lowpass,tuple(FilteredPyramid))
        invtransform = transform.inverse(pyramidd, gain_mask=None)

        IAscombeInv = ((invtransform/2)**2)-(3/8)
        IAscombeInv=numpy.abs(IAscombeInv)
        IAscombeInv[IAscombeInv>numpy.max(Im)]=0
        Images_Filtered.append(IAscombeInv)
    with numpy.errstate(divide='ignore', invalid='ignore'):
        G = numpy.true_divide(Images_Filtered[1],Images_Filtered[0])
        S = numpy.true_divide(Images_Filtered[2],Images_Filtered[0])
    #G=Images_Filtered[1]/Images_Filtered[0]
    #S=Images_Filtered[2]/Images_Filtered[0]
    G=numpy.nan_to_num(G,posinf=0)
    S = numpy.nan_to_num(S,posinf=0)

    original_idxes = numpy.arange(len(G.flatten()))

    return G,S, original_idxes,Images,Images_Filtered

def unmix2species(p3,original_ids,image,P_n,p2):
    """
    Unmixes a mixed FLIM image into two species based on phasor analysis.
    Code used to perform 2 Species STED-FLIM

    Parameters:
    p3 (numpy.ndarray): Phasor coordinates of the mixed FLIM image.
    original_ids (list): List of indices corresponding to the original pixels in the mixed FLIM image.
    image (numpy.ndarray): FLIM image.
    P_n (numpy.ndarray): Phasor centroid coordinates of the first pure species.
    p2 (numpy.ndarray): Phasor centroid coordinates of the second pure species.

    Returns:
    tuple: A tuple containing the unmixed FLIM images for the first species, second species, and the mixing coefficients.
    """
     # Sum all bins of FLIM image to obtain intensity image
    imsum = numpy.sum(image, axis=2)
    #distance between pn and p2
    l= numpy.sum((P_n - p2)**2)
    var=P_n-p2
    #Projection of mixed phasor onto line
    t=numpy.sum((p3 - p2) * (P_n - p2), axis=1) / l
# Clipping the projection between 0 and 1 between P_n and p2
    t = numpy.where(t<0, 0, t)
    t = numpy.where(t>1, 1, t)
    Solve=numpy.array([t,(1-t)])
    # Reassign the values of fractional components to the image
    imsum_flat_lin1 = numpy.array(imsum.flatten()*0, dtype='float32')
    imsum_flat_lin2 = numpy.array(imsum.flatten()*0, dtype='float32')
    for idx,j in enumerate(original_ids) :
        imsum_flat_lin1[j] = t[idx]
        imsum_flat_lin2[j] = 1-t[idx]

    imsum_flat_lin1 = imsum_flat_lin1.reshape(imsum.shape)
    imsum_flat_lin2 = imsum_flat_lin2.reshape(imsum.shape)
    return imsum_flat_lin1,imsum_flat_lin2,Solve


def unmix3species(p3,original_ids,image,P_n,p2,p0):
    """
    Unmixes a mixed phasor into 3 fractional components using a linear system of equations.The fractional components are rescaled so their values are between 0 and 1 between two centroids of the mixed image.
    Code used to perform 2 Species SPLIT-STED

    Args:
        p3 (numpy.ndarray): Mixed phasor array of shape (N, 2), where N is the number of pixels. (g,s coordinates of phasor)
        original_ids (numpy.ndarray): Array of original pixel indices of shape (N,).
        image (numpy.ndarray): FLIM image array of shape (M, N, L), where M is the image dimensions and L is the number of FLIM bins.
        P_n (numpy.ndarray): Control point for the first component of shape (2,).
        p2 (numpy.ndarray): Control point for the second component of shape (2,).
        p0 (numpy.ndarray): Control point for the third component of shape (2,).

    Returns:
        tuple: A tuple containing the unmixing results:
            - imsum_flat_lin1 (numpy.ndarray): Fractional component 1 image array of shape (M, M).
            - imsum_flat_lin2 (numpy.ndarray): Fractional component 2 image array of shape (M, M).
            - imsum_flat_lin3 (numpy.ndarray): Fractional component 3 image array of shape (M, M).
            - Solve (numpy.ndarray): Fractional component values array of shape (3, N).
    """
    
    # Sum all bins of FLIM image to obtain intensity image
    imsum = numpy.sum(image[:,:,10:111], axis=2)
    # Calculate position of 2 centroids of the phasor plot
    kmeans = KMeans(n_clusters = 2, init = 'k-means++', random_state = 42)
    y_kmeans = kmeans.fit_predict(p3)
    CoM_x=kmeans.cluster_centers_[:, 0][:].tolist()
    CoM_y=kmeans.cluster_centers_[:, 1][:].tolist()

    # Project the 2 centroids on the lines connecting the 3 control points
    CoM_x.sort()
    CoM_y.sort(reverse=True) # To get the order CoM = [STED1 < STED2, ...]
    l1 = numpy.sum((P_n - p0) ** 2)  # distance between P_n and p2
    var1 = p0 - P_n
    t_max1 = numpy.sum(([CoM_x[1],CoM_y[1]] - P_n) * (p0 - P_n)) / l1
    l2 = numpy.sum((p2 - p0) ** 2)  # distance between p0 and p2
    var2 = p0 - p2
    t_max2 = numpy.sum(([CoM_x[1],CoM_y[1]] - p2) * (p0 - p2)) / l2
    t_min1=numpy.sum(([CoM_x[0],CoM_y[0]] - P_n) * (p0 - P_n)) / l1
    t_min2=numpy.sum(([CoM_x[0],CoM_y[0]] - p2) * (p0 - p2)) / l2
    pmax1= P_n + numpy.multiply(var1, numpy.transpose(t_max1))
    pmax2 = p2 + numpy.multiply(var2, numpy.transpose(t_max2))
    pmin1= P_n + numpy.multiply(var1, numpy.transpose(t_min1))
    pmin2 = p2 + numpy.multiply(var2, numpy.transpose(t_min2))


    # Create system of equations using the control points
    A = numpy.stack((numpy.array([P_n[0], p2[0], p0[0]]), numpy.array([P_n[1], p2[1], p0[1]]), numpy.ones(3)))
     # Invert the matrix
    Ainv = numpy.linalg.inv(A)
    # Solve the system of equations to obtain the fractional components from the mixed phasor
    p3lin = numpy.hstack((p3, numpy.ones((p3.shape[0], 1))))
    Solve = Ainv.dot(numpy.transpose(p3lin))
 
    # Solve the system of equations to obtain the fractional components corresponding to the 2 centroids of the mixed image
    p3lims=numpy.array([pmin1,pmin2,pmax1,pmax2])
    p3lims=numpy.hstack((p3lims, numpy.ones((4, 1))))
    Solvelims=Ainv.dot(numpy.transpose(p3lims))
    
# Clip the fractional components between 0 and 1
    Solve = numpy.where(Solve<0, 0, Solve)
    Solve = numpy.where(Solve>1, 1, Solve)
    Solvelims = numpy.where(Solvelims<0, 0, Solvelims)
    Solvelims = numpy.where(Solvelims>1, 1, Solvelims)

    Solve=Solve/numpy.sum(Solve,axis=0)
    Solvelims=Solvelims/numpy.sum(Solvelims,axis=0)
    
    #print("Solvelims",Solvelims)

    #smax=numpy.max(Solve[2, :])
    #smin=numpy.min(Solve[2, :])
    #normfact = (smax - smin)
# Rescale the fractional components between the 2 centroids
    smax=numpy.max(Solvelims[2, :])
    smin=numpy.min(Solvelims[2, :])
    normfact = (smax - smin)
    
    
    Solvec=Solve.copy()
    ratio0=Solve[0,:]/(Solve[0,:]+Solve[1,:])
    ratio1=Solve[1,:]/(Solve[0,:]+Solve[1,:])

    Solve[2, :] = (Solve[2, :] - smin) / normfact
    
    #for q,check in enumerate(triangle):
     #   if check==True:
     #       Solve[2,q]=1
    diff=Solve[2,:]-Solvec[2,:]

    Solvec[0,:]=Solve[0,:]-(diff*ratio0)
    Solvec[1,:]=Solve[1,:]-(diff*ratio1)
    Solvec[2,:]=Solve[2,:]
    Solvec[numpy.isnan(Solvec)] = 0
    Solve=Solvec

    Solve = numpy.where(Solve<0, 0, Solve)
    Solve = numpy.where(Solve>1, 1, Solve)
    Solve=Solve/numpy.sum(Solve,axis=0)
    compl_check = numpy.sum(Solve,axis=0)
    print("compl_check",compl_check.min(), compl_check.mean(), compl_check.max())
    print("Solve",Solve.shape, Solve.min(), Solve.max()) 

 # Reassign the values of fractional components to the image
    imsum_flat_lin1 = numpy.array(imsum.flatten() * 0, dtype='float32')
    imsum_flat_lin2 = numpy.array(imsum.flatten() * 0, dtype='float32')
    imsum_flat_lin3 = numpy.array(imsum.flatten() * 0, dtype='float32')

    for idx, j in enumerate(original_ids):
        imsum_flat_lin1[j] = Solve[0, idx]
        imsum_flat_lin2[j] = Solve[1, idx]
        imsum_flat_lin3[j] = Solve[2, idx]
    imsum_flat_lin1 = imsum_flat_lin1.reshape(imsum.shape)
    imsum_flat_lin2 = imsum_flat_lin2.reshape(imsum.shape)
    imsum_flat_lin3 = imsum_flat_lin3.reshape(imsum.shape)
    return imsum_flat_lin1, imsum_flat_lin2,imsum_flat_lin3, Solve

def unmix3species_norescale(p3,original_ids,image,P_n,p2,p0):
    """
    Unmixes a mixed phasor into 3 fractional components using a linear system of equations.

    Args:
        p3 (numpy.ndarray): Mixed phasor array of shape (N, 2), where N is the number of pixels. (g,s coordinates of phasor)
        original_ids (numpy.ndarray): Array of original pixel indices of shape (N,).
        image (numpy.ndarray): FLIM image array of shape (M, N, L), where M is the image dimensions and L is the number of FLIM bins.
        P_n (numpy.ndarray): Control point for the first component of shape (2,).
        p2 (numpy.ndarray): Control point for the second component of shape (2,).
        p0 (numpy.ndarray): Control point for the third component of shape (2,).

    Returns:
        tuple: A tuple containing the unmixing results:
            - imsum_flat_lin1 (numpy.ndarray): Fractional component 1 image array of shape (M, M).
            - imsum_flat_lin2 (numpy.ndarray): Fractional component 2 image array of shape (M, M).
            - imsum_flat_lin3 (numpy.ndarray): Fractional component 3 image array of shape (M, M).
            - Solve (numpy.ndarray): Fractional component values array of shape (3, N).
    """
    # Sum all bins of FLIM image to obtain intensity image
    imsum = numpy.sum(image[:,:,10:111], axis=2)

    # Create system of equations using the control points   
    A = numpy.stack((numpy.array([P_n[0], p2[0], p0[0]]), numpy.array([P_n[1], p2[1], p0[1]]), numpy.ones(3)))

     # Invert the matrix
    Ainv = numpy.linalg.inv(A)

     # Solve the system of equations to obtain the fractional components from the mixed phasor
    p3lin = numpy.hstack((p3, numpy.ones((p3.shape[0], 1))))
    Solve = Ainv.dot(numpy.transpose(p3lin))

    # Clip the fractional components between 0 and 1
    Solve = numpy.where(Solve<0, 0, Solve)
    Solve = numpy.where(Solve>1, 1, Solve)
    Solve=Solve/numpy.sum(Solve,axis=0)

    # Reassign the values of fractional components to the image
    imsum_flat_lin1 = numpy.array(imsum.flatten() * 0, dtype='float32')
    imsum_flat_lin2 = numpy.array(imsum.flatten() * 0, dtype='float32')
    imsum_flat_lin3 = numpy.array(imsum.flatten() * 0, dtype='float32')
    for idx, j in enumerate(original_ids):
        imsum_flat_lin1[j] = Solve[0, idx]
        imsum_flat_lin2[j] = Solve[1, idx]
        imsum_flat_lin3[j] = Solve[2, idx]
    imsum_flat_lin1 = imsum_flat_lin1.reshape(imsum.shape)
    imsum_flat_lin2 = imsum_flat_lin2.reshape(imsum.shape)
    imsum_flat_lin3 = imsum_flat_lin3.reshape(imsum.shape)
    return imsum_flat_lin1, imsum_flat_lin2,imsum_flat_lin3, Solve


def SPLIT_STED(p3, P_n, p2, p3_min, p3_max, image, original_ids):
    """
    Perform SPLIT-STED analysis on a STED-FLIM image.

    Args:
        p3 (numpy.ndarray): Phasor coordinates of STED-FLIM image
        P_n (numpy.ndarray): Reference point corresponding to confocal phasor centroid
        p2 (numpy.ndarray): Second reference point corresponding to short lifetime. most often (1,0)
        p3_min (numpy.ndarray): phasor coordinate of first centroid of STED-FLIM image phasor.
        p3_max (numpy.ndarray): phasor coordinate of first centroid of STED-FLIM image phasor.
        image (numpy.ndarray): Input image.
        original_ids (numpy.ndarray): Original indices of p3 phasor points in image space

    Returns:
        tuple: A tuple containing:
            - im_fract (numpy.ndarray): Image containing fractional components for each pixel
            - projection (numpy.ndarray): Phasor coordinates of the projected points.
            - t2 (numpy.ndarray): Fractional component values

    """
    l2 = numpy.sum((P_n - p2) ** 2)  # distance between P_n and p2
    var = p2 - P_n

    #Project phasor on line connecting p2 and P_n (Parametric coordinates)
    t = numpy.sum((p3 - P_n) * (p2 - P_n), axis=1) / l2 
#Project phasor centroids on line connecting p2 and P_n
    t_min = numpy.sum((p3_min - P_n) * (p2 - P_n)) / l2  
    t_max = numpy.sum((p3_max - P_n) * (p2 - P_n)) / l2
    normfact = (t_max - t_min)
    #Rescale projected phasor points between 2 centroids so fractions will be [0,1]
    t2 = numpy.array([t, t])
    t2 = (t2 - t_min) / normfact
    t2 = numpy.where(t2 < 0, 0, t2)
    t2 = numpy.where(t2 > 1, 1, t2)

# Find g,s coordinates of projected points (for phasor space graphs)
    projection = P_n + numpy.multiply(p2 - P_n, numpy.transpose(t2)) 
 # reassign fractions to pixels 
    im_fract = numpy.array(numpy.sum(image,axis=2).flatten()*0, dtype='float32')
    for idx, j in enumerate(original_ids):
        im_fract[j] = 1-t2[0, idx]
    im_fract = im_fract.reshape(image.shape[0:2])
    return im_fract,projection,t2
