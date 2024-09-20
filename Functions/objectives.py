
"""This module contains classes that implement several objectives to optimize.
One can define a new objective by inheriting abstract class :class:`Objective`.
"""

from abc import ABC, abstractmethod

import numpy
import itertools
import warnings

from statsmodels.tsa.stattools import acf
from scipy.ndimage import gaussian_filter
from scipy import optimize
from skimage.transform import resize
from sklearn.metrics import mean_squared_error

#from . import fsc
# import src.utils as utils
# import src.user as user

import decorr
# from . import utils
# from . import user

from scipy.ndimage import gaussian_filter1d


class Objective(ABC):
    """Abstract class to implement an objective to optimize. When inheriting this class,
    one needs to define an attribute `label` to be used for figure labels, and a
    function :func:`evaluate` to be called during optimization.
    """
    @abstractmethod
    def evaluate(self, sted_stack, confocal_init, confocal_end, sted_fg, confocal_fg):
        """Compute the value of the objective given the result of an acquisition.

        :param sted_stack: A list of STED images.
        :param confocal_init: A confocal image acquired before the STED stack.
        :param concofal_end: A confocal image acquired after the STED stack.
        :param sted_fg: A background mask of the first STED image in the stack
                        (2d array of bool: True on foreground, False on background).
        :param confocal_fg: A background mask of the initial confocal image
                            (2d array of bool: True on foreground, False on background).
        """
        raise NotImplementedError

    def mirror_ticks(self, ticks):
        """Tick values to override the true *tick* values for easier plot understanding.

        :param ticks: Ticks to replace.

        :returns: New ticks or None to keep the same.
        """
        return None


#class Signal_Ratio(Objective):
#    """Objective corresponding to the signal to noise ratio (SNR) defined by
#
#    .. math::
#        \\text{SNR} = \\frac{\\text{STED}_{\\text{fg}}^{75} - \overline{\\text{STED}_{\\text{fg}}}}{\\text{Confocal1}_{\\text{fg}}^{75}}
#
#    where :math:`\\text{image}^q` and :math:`\\overline{\\text{image}}` respectively
#    denote the :math:`q`-th percentile signal on an image and the average signal
#    on an image, and :math:`\\text{STED}_{\\text{fg}}`, :math:`\\text{Confocal1}_{\\text{fg}}`, and
#    :math:`\\text{Confocal2}_{\\text{fg}}` respectively refer to the foreground of the STED image
#    and confocal images acquired before and after.
#
#    :param float percentile: :math:`q`-th percentile in :math:`[0,100]`.
#    """
#    def __init__(self, percentile):
#        self.label = "Signal Ratio"
#        self.select_optimal = numpy.argmax
#        self.percentile = percentile
#
#    def evaluate(self, sted_stack, confocal_init, confocal_end, sted_fg, confocal_fg):
#        """Compute the signal to noise ratio (SNR) given the result of an acquisition.
#
#        :param sted_stack: A list of STED images.
#        :param confocal_init: A confocal image acquired before the STED stack.
#        :param concofal_end: A confocal image acquired after the STED stack.
#        :param sted_fg: A background mask of the first STED image in the stack
#                        (2d array of bool: True on foreground, False on background).
#        :param confocal_fg: A background mask of the initial confocal image
#                            (2d array of bool: True on foreground, False on background).
#
#        :returns: :math:`0` if no STED foreground, None if :math:`\\text{SNR} < 0` (error), or
#                  SNR value otherwise.
#
#        """
#        if numpy.any(sted_fg):
#            foreground = numpy.percentile(sted_stack[0][sted_fg], self.percentile)
#            background = numpy.mean(sted_stack[0][numpy.invert(sted_fg)])
#            ratio = (foreground - background) / numpy.percentile(confocal_init[confocal_fg], self.percentile)
#            if ratio < 0:
#                return None
#            else:
#                return ratio
#        else:
#            return 0


class Signal_Ratio(Objective):
    def __init__(self, norm=True):
        self.label = "Signal Ratio"
        self.select_optimal = numpy.argmax
        self.norm = norm

    def evaluate(self, sted_stack, confocal_init, confocal_end, sted_fg, confocal_fg):
        if numpy.any(sted_fg):
            sted = sted_stack[0]
            I = numpy.median(sted[sted_fg]) # I=S+B
            B = numpy.mean(sted[numpy.invert(sted_fg)])
            SNR = (I-B)/numpy.sqrt(I)
            if self.norm:
                I = numpy.median(confocal_init[confocal_fg]) # I=S+B
                B = numpy.mean(confocal_init[numpy.invert(confocal_fg)])
                SNR =  SNR / ( (I-B)/numpy.sqrt(I))
            return SNR
        else:
            return 0
            

class FWHM(Objective):
    """Objective corresponding to the full width at half maximum (FWHM) defined by

    .. math::
        \\text{average}(|2.3558 \\cdot \\sigma |) \\cdot p \\cdot 10^9

    where :math:`p` is the size of a pixel (in nm) in the STED image, and :math:`\\sigma`
    is the standard deviation estimated from a Gaussian fit (see function
    :func:`utils.gaussian_fit`) on at least three line profiles using function
    :func:`user.get_lines`.

    :param pixelsize: Size of a pixel in a STED image (in nm).
    :param `**kwargs`: This method also takes the keyword arguments for :func:`user.get_lines`.
    """
    def __init__(self, pixelsize, **kwargs):
        self.label = "FWHM (nm)"
        self.select_optimal = numpy.argmin
        self.pixelsize = pixelsize
        self.kwargs = kwargs

    def evaluate(self, sted_stack, confocal_init, confocal_end, sted_fg, confocal_fg):
        """Compute the full width at half maximum (FWHM) given the result of an acquisition.
        It relies on the function :func:`user.get_lines` to request the user to select line
        profiles in the first STED image of the stack. If the user does not select any lines
        in the STED image, ask the user to select line profiles in the initial confocal.

        :param sted_stack: A list of STED images.
        :param confocal_init: A confocal image acquired before the STED stack.
        :param concofal_end: A confocal image acquired after the STED stack.
        :param sted_fg: A background mask of the first STED image in the stack
                        (2d array of bool: True on foreground, False on background).
        :param confocal_fg: A background mask of the initial confocal image
                            (2d array of bool: True on foreground, False on background).

        :returns: The averaged FWHM (in nm) if success, else None.
        """
        lines = user.get_lines(sted_stack[0], 3, minlen=4, deltas=[-1, 0, 1], **self.kwargs)
        if not lines:
            lines = user.get_lines(confocal_init, 3, minlen=4, deltas=[-1, 0, 1], **self.kwargs)
        fwhms = []
        for positions, profile in lines:
            popt = utils.gaussian_fit(positions, profile)
            if popt is not None:
                fwhms.append(numpy.abs(2.3548 * popt[-1]))
        if fwhms: return numpy.mean(fwhms)*self.pixelsize*1e9
        else: return None



class Autocorrelation(Objective):
    """Objective corresponding to the autocorrelation defined as the difference between
    the value at the first maximum and the value at the first minimum following the first
    maximum, for

    :param `**kwargs`: This method also takes the keyword arguments for :func:`user.get_lines`.
    """
    def __init__(self, **kwargs):
        self.label = "Autocorrelation"
        self.select_optimal = numpy.argmax
        self.kwargs = kwargs

    def evaluate(self, sted_stack, confocal_init, confocal_end, sted_fg, confocal_fg):
        lines = user.get_lines(sted_stack[0], minlen=40, deltas=[-1, 0, 1], **self.kwargs)
        profiles = [l[1] for l in lines]
        autocorr = acf(profiles)
        min_val, min_idx = utils.find_first_min(autocorr)
        max_val, max_idx = utils.find_first_max(autocorr, min_idx)
        assert max_val >= min_val
        if max_idx < min_idx:
            return None
        else:
            return max_val - min_val


class Score(Objective):
    """Objective corresponding to the autocorrelation defined as the difference between
    the value at the first maximum and the value at the first minimum following the first
    maximum, for

    :param `**kwargs`: This method also takes the keyword arguments for :func:`user.give_score`.
    """
    def __init__(self, label, select_optimal=numpy.argmax, idx=0, **kwargs):
        self.label = label
        self.select_optimal = select_optimal
        self.idx = idx
        self.kwargs = kwargs

    def evaluate(self, sted_stack, confocal_init, confocal_end, sted_fg, confocal_fg):
        return user.give_score(confocal_init, sted_stack[self.idx], self.label, **self.kwargs)


class Bleach(Objective):
    def __init__(self):
        self.label = "Bleach"
        self.select_optimal = numpy.argmin

    def evaluate(self, sted_stack, confocal_init, confocal_end, sted_fg, confocal_fg):
        signal_i = numpy.mean(confocal_init[confocal_fg])
        signal_e = numpy.mean(confocal_end[confocal_fg])
        bleach = (signal_i - signal_e) / signal_i
        return bleach
        
class Bleach_simulation(Objective):
    def __init__(self):
        self.label = "Bleach simulation"
        self.select_optimal = numpy.argmin

    def evaluate(self, datamap, bleached_datamap, sted, ):
        signal_i = numpy.mean(datamap[datamap!=0])
        signal_e = numpy.mean(bleached_datamap[datamap!=0])
        bleach = (signal_i - signal_e) / signal_i
        return bleach


class ScoreNet(Objective):
    def __init__(self, label, net, select_optimal=numpy.argmax, idx=0):
        self.label = label
        self.net = net
        self.select_optimal = select_optimal
        self.idx = idx

    def evaluate(self, sted_stack, confocal_init, confocal_end, sted_fg, confocal_fg):
        score = self.net.predict(utils.img2float(sted_stack[self.idx]))
        print("Net", self.label, "score", score)
        return score


class FRC(Objective):
    def __init__(self, pixelsize, res_cap=350):
        self.label = "FRC"
        self.select_optimal = numpy.argmax
        self.pixelsize = pixelsize # µm
        self.max_spatialfreq = 1 / (2 * pixelsize) # 1/µm
        self.res_cap=res_cap

    def evaluate(self, sted_stack, confocal_init, confocal_end, sted_fg, confocal_fg):
        sted = numpy.array(sted_stack[0])
        # verify that the STED image is of squared shape
        assert sted.shape[0] == sted.shape[1],\
            "The STED image is not a square, you cannot evaluate the Fourier Ring Correlation!"
        if len(sted_stack)==1:
            imgs = fsc.split_image_array(sted, 2)
        else:
            imgs = [sted, numpy.array(sted_stack[1])]
        fourierringcorr, sigmacurves, freqs = [], [], []
        for im1, im2 in itertools.combinations(imgs, 2):
            frc, nPx = fsc.fourier_shell_corr(im1, im2)
            sigma = fsc.sigma_curve(nPx)
            freq = numpy.arange(frc.shape[0]) / (im1.shape[0] * self.pixelsize)
            fourierringcorr.append(frc)
            sigmacurves.append(sigma)
            freqs.append(freq)

        frc = numpy.mean(numpy.array(fourierringcorr), axis = 0)
        sigma = numpy.mean(numpy.array(sigmacurves), axis = 0)
        freq = numpy.mean(numpy.array(freqs), axis = 0)

        spatialfreq = fsc.meeting_point(fsc.moving_average(frc, 3), freq, fsc.moving_average(sigma, 3))
        if spatialfreq!=0:
            return min(self.res_cap, (1/spatialfreq)/1e-9) #/ self.max_spatialfreq
        else:
            return self.res_cap

    def mirror_ticks(self, ticks):
        return ["{:0.0f}".format(1e+3 / (self.max_spatialfreq * x)) if x > 0 else "" for x in ticks]


class Resolution(Objective):
    def __init__(self, pixelsize, res_cap=350):
        self.label = "Resolution (nm)"
        self.select_optimal = numpy.argmin
        self.pixelsize = pixelsize
#            self.kwargs = kwargs
        self.res_cap=res_cap

    def evaluate(self, sted_stack, confocal_init, confocal_end, sted_fg, confocal_fg):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = decorr.calculate(sted_stack[0])*self.pixelsize/1e-9
        if res > self.res_cap:
            res = self.res_cap
        return res

class Squirrel(Objective):
    """
    Implements the `Squirrel` objective

    :param method: A `str` of the method used to optimize
    :param normalize: A `bool` wheter to normalize the images
    """
    def __init__(self, method="L-BFGS-B", normalize=False):

        self.method = method
        self.bounds = (-numpy.inf, numpy.inf), (-numpy.inf, numpy.inf), (0, numpy.inf)
        self.x0 = (1, 0, 1)
        self.normalize = normalize
        self.select_optimal = numpy.argmin


    def return_map(self, sted_stack, confocal_init, confocal_end, sted_fg, confocal_fg):
        result = self.optimize(sted_stack[0], confocal_init)

        alpha, beta, sigma = result.x
        super_resolution, reference = sted_stack[0], confocal_init
        convolved = self.convolve(super_resolution, alpha, beta, sigma)
        if self.normalize:
            reference = (reference - reference.min()) / (reference.max() - reference.min() + 1e-9)
            convolved = (convolved - convolved.min()) / (convolved.max() - convolved.min() + 1e-9)
    
        errmap=(convolved -reference) ** 2
        return errmap, result.x
       

    def evaluate(self, sted_stack, confocal_init, confocal_end, sted_fg, confocal_fg):
        """
        Evaluates the objective

        :param sted_stack: A list of STED images.
        :param confocal_init: A confocal image acquired before the STED stack.
        :param concofal_end: A confocal image acquired after the STED stack.
        :param sted_fg: A background mask of the first STED image in the stack
                        (2d array of bool: True on foreground, False on background).
        :param confocal_fg: A background mask of the initial confocal image
                            (2d array of bool: True on foreground, False on background).
        """
        # Optimize
        if not numpy.any(sted_stack[0]):
            return mean_squared_error(confocal_init[confocal_fg], sted_stack[0][confocal_fg], squared=True)
        
        # Optimize

        result = self.optimize(sted_stack[0], confocal_init)


        return self.squirrel(result.x, sted_stack[0], confocal_init)

    def squirrel(self, x, *args):
        """
        Computes the reconstruction error between
        """
        alpha, beta, sigma = x

        super_resolution, reference = args
        convolved = self.convolve(super_resolution, alpha, beta, sigma)
        if self.normalize:
            reference = (reference - reference.min()) / (reference.max() - reference.min() + 1e-9)
            convolved = (convolved - convolved.min()) / (convolved.max() - convolved.min() + 1e-9)
        
        error = mean_squared_error(reference, convolved, squared=True)
        

        return error


    def optimize(self, super_resolution, reference):
        """
        Optimizes the SQUIRREL parameters

        :param super_resolution: A `numpy.ndarray` of the super-resolution image
        :param reference: A `numpy.ndarray` of the reference image

        :returns : An `OptimizedResult`
        """


        result = optimize.minimize(
                self.squirrel, self.x0, args=(super_resolution, reference),
                method="L-BFGS-B", bounds=((-numpy.inf, numpy.inf), (-numpy.inf, numpy.inf), (0, numpy.inf)))

        return result

    def convolve(self, img, alpha, beta, sigma):
        """
        Convolves an image with the given parameters
        """
        return gaussian_filter(img * alpha + beta, sigma=sigma)

