
"""This module contains classes that implement several objectives to optimize.
One can define a new objective by inheriting abstract class :class:`Objective`.
"""

from abc import ABC, abstractmethod

import numpy
import itertools
import warnings
from scipy.ndimage import gaussian_filter
from scipy import optimize
from sklearn.metrics import mean_squared_error
import decorr



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


class Bleach(Objective):
    def __init__(self):
        self.label = "Bleach"
        self.select_optimal = numpy.argmin

    def evaluate(self, sted_stack, confocal_init, confocal_end, sted_fg, confocal_fg):
        signal_i = numpy.mean(confocal_init[confocal_fg])
        signal_e = numpy.mean(confocal_end[confocal_fg])
        bleach = (signal_i - signal_e) / signal_i
        return bleach
        


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

