
import numpy
import tifffile
import matplotlib
import numpy
from scipy import ndimage
from matplotlib import pyplot, colors

def mean_filter(array, footprint=(3,3)):
    """
    Smooths an `numpy.ndarray` by applying a mean filter

    :param ary: A `numpy.ndarray`
    :param footprint: A `tuple` of the mean filter size
    """
    kernel = numpy.ones(footprint)
    return ndimage.convolve(array, kernel) / kernel.sum()

class LifetimeOverlayer:
    """
    Creates a `LifetimeOverlayer`. This allows a user to map a lifetime image
    with an intensity image. At each position, the color is weighted by its
    intensity in an intensity image.
    """
    def __init__(self, lifetime, intensity=None, cname="rainbow"):
        """
        Instantiates the `LifetimeOverlayer`

        :param lifetime: A 2D `numpy.ndarray` of the lifetime
        :param intensity: A 2D/3D `numpy.ndarray` of the intensity image
        :param cname: A `str` of the colormap name
        """
        self.lifetime = lifetime
        if isinstance(intensity, type(None)):
            self.intensity = numpy.ones_like(self.lifetime)
        else:
            self.intensity=intensity
        #     self.intensity = self.verify_intensity(intensity)
        self.cname = cname

    def get_overlay(self, lifetime_minmax=(0., 5.), intensity_minmax=(0., 1.)):
        """
        Computes the lifetime overlay

        :param lifetime_minmax: A `tuple` of the lifetime minimum and maximum
        :param intensity_minmax: A `tuple` of the intensity minimum and maxium

        :returns : A `numpy.ndarray` of the weighted colors
        """
        # Normalize intensity image
        _min, _max = intensity_minmax
        intensity = numpy.clip(
            (self.intensity - _min) / (_max - _min), 0, 1
        )
    
        # Normalize lifetime image
        _min, _max = lifetime_minmax
        norm = colors.Normalize(vmin=_min, vmax=_max)
        cmap = matplotlib.cm.ScalarMappable(norm, self.cname)
        lifetime_rgb = cmap.to_rgba(self.lifetime)[:, :, :3]

        # Convert to hsv and applies intensity mapping
        lifetime_hsv = colors.rgb_to_hsv(lifetime_rgb)
        lifetime_hsv[:, :, -1] = intensity

        lifetime_rgb = colors.hsv_to_rgb(lifetime_hsv)
        return lifetime_rgb, cmap
    
    def get_overlay_test(self, lifetime_minmax=(0., 5.), intensity_minmax=(0., 1.)):
        """
        Computes the lifetime overlay

        :param lifetime_minmax: A `tuple` of the lifetime minimum and maximum
        :param intensity_minmax: A `tuple` of the intensity minimum and maxium

        :returns : A `numpy.ndarray` of the weighted colors
        """
        # Normalize intensity image
        _min, _max = intensity_minmax
        intensity = numpy.clip(
            (self.intensity - _min) / (_max - _min), 0, 1
        )
    
        # Normalize lifetime image
        _min, _max = lifetime_minmax
        norm = colors.Normalize(vmin=_min, vmax=_max)
        #cmap = matplotlib.cm.ScalarMappable(norm, self.cname)

        red=self.lifetime[:,:,1]+self.lifetime[:,:,2]
        green=self.lifetime[:,:,0]+self.lifetime[:,:,2]
        blue=self.lifetime[:,:,1]+self.lifetime[:,:,0]
        
        lifetime_rgb =numpy.dstack((red,green,blue))
        lifetime_rgb=numpy.where(lifetime_rgb>1, 1,lifetime_rgb)
        print(lifetime_rgb.shape)
        #lifetime_rgb = cmap.to_rgba(self.lifetime)[:, :, :3]

        # Convert to hsv and applies intensity mapping
        lifetime_hsv = colors.rgb_to_hsv(lifetime_rgb)
        lifetime_hsv[:, :, -1] = intensity

        lifetime_rgb = colors.hsv_to_rgb(lifetime_hsv)
        return lifetime_rgb

    def verify_intensity(self, intensity):
        """
        Ensures that the intensity map is 2D and minimum value is 0.

        :param intensity: A `numpy.ndarray` of the intensity

        :returns : A `numpy.ndarray` of the intensity map
        """
        if intensity.min() == 0:
            print(0)
            return intensity
        if intensity.dtype == numpy.float8:
             print('float8')
             return intensity - 2**8 / 2
        elif intensity.dtype == numpy.float16:
            print('float16')
            return intensity - 2**16 / 2
        elif intensity.dtype == numpy.float32:
            print('float32')
            return intensity - 2**32 / 2
        else:
            print(intensity.dtype)
            return intensity

        
        # if intensity.min() != 0:
        #     # Assumes intensity is uint16
        #     print(intensity.min(), intensity.max())
        #     intensity = intensity - 2 ** 15
        #     print(intensity.min(), intensity.max())
        # if intensity.ndim > 2:
        #     # Assumes intensity image is (C, H, W)
        #     intensity = intensity.sum(axis=0)
        # return intensity

'''Exemple
if __name__ == "__main__":

    lifetime = tifffile.imread("lt_map.tif")
    lifetime = mean_filter(lifetime)

    intensity = tifffile.imread("star635_lifetime_golabl.tif")

    overlayer = LifetimeOverlayer(lifetime, intensity)
    lifetime_rgb, cmap = overlayer.get_overlay(
        lifetime_minmax=(1.0, 4.5),
        intensity_minmax=(0., 1.)
    )

    fig, ax = pyplot.subplots()
    ax.imshow(lifetime_rgb)
    cbar = pyplot.colorbar(cmap, ax=ax)
    cbar.set_label("Lifetime (ns)")
    pyplot.show()
'''