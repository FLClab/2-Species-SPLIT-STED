# 2 Species SPLIT-STED

All the codes used to produce figures and analyze data for:  
Enhancing STED Microscopy via Fluorescence Lifetime Unmixing and Filtering in Two-Species SPLIT-STED.

The Confocal- and STED-FLIM images of neuronal proteins dataset is available to downlad from the [Zenodo dataset](https://doi.org/10.5281/zenodo.15438494)
## Installation and Requirements

The source code `2 Species SPLIT-STED` relies on Python scientific librairies. The source code was tested in a Python 3.11 environnement. We provide a `requirements.txt` file to facilitate the installation of the necessary dependencies.

Assuming the users have a working installation of Python on their computer (we recommend using [Anaconda](https://docs.anaconda.com/anaconda/install/)), the users should create a new Python 3.11 environnement to avoid impacting on other file dependencies. 

```bash
conda create -n FLIM python=3.11.4
conda activate FLIM
pip install -r requirements.txt
```
If available, install the specpy package provided with the Imspector software for your specific python version. All data is provided as both .tiff and .msr files and the scripts can read both. If you want to use these scripts with other data formats, simply change the load_image and select_channel functions.

## Folder contents
### Acquisition
Codes to perform automatic image acquisition with different depletion powers on an Abberior Expert Line STED microscope. Based on the [specpy](https://pypi.org/project/specpy/) and [Abberior-STED](https://https://github.com/FLClab/Abberior-STED) libraries

- `FLIM_AutoAcquire_VaryPower_ConfocalsPrePost.py`: Main program that coordinates the acquisition sequence, randomly selects depletion power values from a list of repeated values and sets the value for acquisition. For each region, it lauches the acquisition a predefined sequence of images :
1) Confocal image
2) Pair of STED-FLIM and Confocal-FLIM images
3) Confocal image.
### Functions
Functions called in other folders
- `Main_Functions.py`: Contains phasor calibration codes, foreground detection
- `Phasor_functions.py`: Contains code to calculate both median and [CWF](https://doi.org/10.1364/BOE.420953) filtered phasor distributions, codes to perform unmixing and SPLIT-STED. 
- `LineProfile_Nchannels_tifffile_1graphperline_.py`: Script to plot intensity profiles (with interactive window to select lines)
- `objectives.py` : Contains photobleaching and [SQUIRREL metric](https://doi.org/10.1038/nmeth.4605) calculation functions
- `decorr.py`: Function to calculate image resolution based on [Decorrelation analysis](https://doi.org/10.1038/s41592-019-0515-7). Script originally made for [this paper](https://doi.org/10.1101/2024.03.25.586697) and available [here](https://github.com/FLClab/bandit-optimization)


###  Histograms
Curve fitting of FLIM histograms
- `IRF measurement_Gaussfit.py`: Fits a gaussian function onto the histogram of an IRF measruement (imaging a sample of gold nanoparticles)
- `Hist_monoexp_MLE_foreground.py`: Sums histograms of all pixels in the foreground of a FLIM image and fits a single exponential function to the data using MLE error estimation
- `Hist_monoexp_MLE_foreground_AllFolder.py`: Same as **Hist_monoexp_MLE_foreground** but repeats for all the FLIM images in a folder and saves results to a csv file.
- `Hist_monoexp_MLE_pixelwise.py`: Fits a single exponential function to the histogram of every pixel in an image data using MLE error estimation.
- `Hist_monoexp_MLE_pixelwise_AllFolder.py`: Same as **Hist_monoexp_MLE_pixelwise** but repeats for all the FLIM images in a folder

### PhasorDistribution
Display and analysis of phasor distributions of STED-FLIM images

- `PhasorDistribution_STED_EllipseFit_2Species.py` Analysis of phasor distribution trajectories as a function of STED depletion power by fitting ellipses on phasor distributions. Saves metrics (shortest distance between ellipses, IOU, distance between ellipse centroids , ...) to a csv file and graphs properties (mean phasor centroid, ellipse dimensions) for each depletion power.
- `ReadCSVtoGraphs_EllipseFit.py`: Reads csv file produced by **PhasorDistribution_STED_EllipseFit_2Species** and plots the different metrics as a function of depletion power
- `PhasorDistribution_STED_3DGraph.py` : Makes 3D graph of phasor distributions of 2 stainings as a function of STED depletion power 
- `PhasorDistribution_PrevsPostCalibration.py`: Plots a phasor before and after the calibration with the IRF
- `PhasorDistribution_RawvsMedianvsCWF.py`: Plots the same phasor with different filtering techniques (raw,median and CWF)

### SPLIT-STED
1 Species [SPLIT-STED](https://doi.org/10.1039/C8NR07485B) to improve resolution of single fluorophore STED-FLIM images
- `SPLIT-STED_CWF_multiSTEDpercent_Metrics.py`: 1 species SPLIT-STED of all images in a folder with [*CWF filtering*](https://doi.org/10.1364/BOE.420953) of phasors. Saves resolution, photobleaching and SQUIRREL metrics of the input and output images in a csv file
- `SPLIT-STED_Median_multiSTEDpercent_Metrics.py`: 1 species SPLIT-STED of all images in a folder with *median filtering* of phasors. Saves resolution, photobleaching and SQUIRREL metrics of the input and output images in a csv file
- `ReadCSVstoGraphs_EllipsesvsSPLITresolution.py`: Read csv files generated by SPLIT-STED codes and `PhasorDistribution_STED_EllipseFit_2Species.py` and produces graphs


### Unmixing
Algorithms to separate dyes of different lifetimes in Confocal-FLIM and STED-FLIM images using [linear systems of equations in phasor space](https://doi.org/10.1088/2050-6120/ab8570)
- `Unmixing_2SpeciesSPLITSTED.ipynb`is a Jupyter Notebook implementation of 2 Species SPLIT-STED. A small set of example data is downloaded by the script from this repository into an **Example_data** subfolder if it is not already present.
- `Unmixing_2SpeciesConfocalFLIM__allfolder.py` : Separate phasors of Confocal-FLIM images of double-stained samples into two fractional components. Uses the same pair of control images for the entire folder of mixed images
- `Unmixing_2SpeciesSTEDFLIM_allfolder.py`: Separate phasors of STED-FLIM images of double-stained samples into two fractional components. Uses the same control images for all mixture images acquired with the same depletion power (1 pair of controls per depletion power)
- `Unmixing_2SpeciesSPLITSTED_allfolder.py` :  Separate phasors of STED-FLIM images of double-stained samples into three fractional components. Uses the same control images (4 per dye, confocal and 3 depletion powers) for all mixture images.

### Simulation
Creates simulated double-staining images by summing single-staining STED-FLIM images and performs different unmixing algorithms to evaluate their performance in a controlled setting.
- `HistogramMLE_MonoExponential_Simulationimages.py`: Generate simulated images and perform MLE histogram fitting to generate color-coded lifetime images of the ground truth and mixed images
- `Simulation_2SpeciesSTED-FLIM_looppowers.py` : Generate simulated images and perform 2 Species STED-FLIM unmixing. Calculate metrics (resolution, squirrel) and save them to a csv file. Also save input images and unmixing results (color-coded phasor and images) for each simulated combination.
- `Simulation_2SpeciesSPLIT-STED_looppowers.py`: Generate simulated images and perform 2 Species SPLIT-STED unmixing. Calculate metrics (resolution, squirrel) and save input images and unmixing results (color-coded phasor and images)
- `ReadCSVstoGraphs_Simulation.py` : Read the csv files produced by **Simulation_2SpeciesSTED-FLIM_looppowers.py` and **Simulation_2SpeciesSPLIT-STED_looppowers.py` and produce comparison graphs of the performance metrics of both unmixing methods as a function of the depletion power.

