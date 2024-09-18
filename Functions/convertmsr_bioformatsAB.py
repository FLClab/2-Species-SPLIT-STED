"""
If running into import errors:
- Remove the try catch block to see the error --> probably "JVM not found"
- Make sure you have a working java version installed:
    e.g.:
        java version "1.8.0_311"
        Java(TM) SE Runtime Environment (build 1.8.0_311-b11)
        Java HotSpot(TM) 64-Bit Server VM (build 25.311-b11, mixed mode)
- Make sure you have numpy installed (recommended to work in a virtual env):
    $ pip install numpy
- Uninstall and reinstall packages, starting with javabridge:
    $ pip uninstall javabridge
    $ pip uninstall python-bioformats
    $ pip install javabridge
    $ pip install python-bioformats
"""


import os
import numpy

try:
    import javabridge
    import bioformats
except ImportError:
    print("Bioformats does not seem to be installed on your machine...")
    print("Try running `pip install python-bioformats`")
    exit()

class MSRReader:
    """
    Creates a `MSRReader`. It will take some time to create the object

    :param logging_level: A `str` of the logging level to use {WARN, ERROR, OFF}

    :usage :
        with MSRReader() as msrreader:
            data = msrreader.read(file)
            image = data["STED_640"]
    """
    def __init__(self, logging_level="OFF"):

        # Starts the java virtual machine
        javabridge.start_vm(class_path=bioformats.JARS)

        rootLoggerName = javabridge.get_static_field("org/slf4j/Logger","ROOT_LOGGER_NAME", "Ljava/lang/String;")
        rootLogger = javabridge.static_call("org/slf4j/LoggerFactory","getLogger", "(Ljava/lang/String;)Lorg/slf4j/Logger;", rootLoggerName)
        logLevel = javabridge.get_static_field("ch/qos/logback/classic/Level", logging_level, "Lch/qos/logback/classic/Level;")
        javabridge.call(rootLogger, "setLevel", "(Lch/qos/logback/classic/Level;)V", logLevel)

    def read(self, msrfile):
        """
        Method that implements a `read` of the given `msrfile`

        :param msrfile: A file path to a `.msr` file

        :returns : A `dict` where each keys corresponds to a specific image
                   in the measurement file
        """
        data = {}
        with bioformats.ImageReader(path=msrfile) as reader:
            metadata = bioformats.OMEXML(bioformats.get_omexml_metadata(path=reader.path))

            # Retreives the number of series
            series = metadata.get_image_count()

            # We iterate over each serie
            rdr = reader.rdr
            for serie in range(series):
                rdr.setSeries(serie)
                X, Y, Z, T, C = rdr.getSizeX(), rdr.getSizeY(), rdr.getSizeZ(), rdr.getSizeT(), rdr.getSizeC()
                Zs = []
                for z in range(Z):
                    Ts = []
                    for t in range(T):
                        Cs = []
                        for c in range(C):
                            image = reader.read(z=z, t=t, c=c, series=serie, rescale=False)
                            Cs.append(image)
                        Ts.append(Cs)
                    Zs.append(Ts)

                # Avoids single axes in data
                image = numpy.array(Zs).squeeze()

                # Stores in data folder
                image_metadata = metadata.image(serie)
                data[image_metadata.get_Name()] = image
        return data

    def close(self):
        javabridge.kill_vm()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __del__(self):
        self.close()

if __name__ == "__main__":

    msrfiles = [
        '/Users/marielafontaine/valeria-s3/flclab-abberior-sted/mlafontaine/2022-05-17/1_color_Cy3_template_mCardinal.msr',
        '/Users/marielafontaine/valeria-s3/flclab-abberior-sted/mlafontaine/2022-05-17/Tetraspeck_Alignment_.msr', 
        '/Users/marielafontaine/valeria-s3/flclab-abberior-sted/mlafontaine/2022-05-17/Tetraspeck_Alignment_FLIM.msr'
    ]
# Warning : cette etape prend du temps
    with MSRReader() as msrreader:
        for msrfile in msrfiles:
            data = msrreader.read(msrfile)
            for key, value in data.items():
                print(key, value.shape)
