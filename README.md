# 2 Species SPLIT-STED

All the codes used to produce figures and analyze data for:  
Three-Component Phasor-Based Linear Unmixing Approach for Simultaneous Fluorophore Discrimination and Resolution Improvement of Stimulated Emission Depletion Microscopy Images Using Fluorescence Lifetime.

Dataset is available to downlad from the [website](https://s3.valeria.science/flclab-2-species-split-sted/index.html) 
## Installation and Requirements

The source code `2 Species SPLIT-STED` relies on Python scientific librairies. The source code was tested in a Python 3.11 environnement. We provide a `requirements.txt` file to facilitate the installation of the necessary dependencies.

Assuming the users have a working installation of Python on their computer (we recommend using [Anaconda](https://docs.anaconda.com/anaconda/install/)), the users should create a new Python 3.11 environnement to avoid impacting on other file dependencies. 

```bash
conda create -n FLIM python=3.11.4
conda activate FLIM
pip install -r requirements.txt
```
If available, install the specpy package provided with the Imspector software for your specific python version. All data is provided as tiff images and the scripts can read both .tiff and msr files. If you want to use these scripts with other data formats, simply change the load_image and select_channel functions.

## Folder contents
### Acquisition
Codes to perform automatic image acquisition with different depletion powers on an Abberior Expert Line STED microscope. Based on the [specpy](https://pypi.org/project/specpy/) and [Abberior-STED](https://https://github.com/FLClab/Abberior-STED) libraries

* **FLIM_AutoAcquire_VaryPower_ConfocalsPrePost**: Main program that coordinates the acquisition sequence, randomly selects depletion power values from a list of repeated values and sets the value for acquisition. For each region, it lauches the acquisition a predefined sequence of images :
1) Confocal image
2) Pair of STED-FLIM and Confocal-FLIM images
3) Confocal image.
### Functions
Functions called in other folders
* **Main_Functions**: Contains phasor calibration codes, foreground detection

<img src="images/Calibration_Phasor.png" width="320"/>

* **Phasor_functions**: Contains code to calculate both median and [CWF](https://doi.org/10.1364/BOE.420953) filtered phasor distributions, codes to perform unmixing and SPLIT-STED. 
* **LineProfile_Nchannels_tifffile_1graphperline_**: Script to plot intensity profiles (with interactive window to select lines)
* **objectives** : Contains photobleaching and [SQUIRREL metric](https://doi.org/10.1038/nmeth.4605) calculation functions
* **decorr** : Function to calculate image resolution based on [Decorrelation analysis](https://doi.org/10.1038/s41592-019-0515-7). Script originally made for [this paper](https://doi.org/10.1101/2024.03.25.586697) and available [here](https://github.com/FLClab/bandit-optimization)
* **lifetime** : Function to create images colorcoded by lifetime overlayed on the intensity values.

###  Histograms
Curve fitting of FLIM histograms
* **IRF measurement_Gaussfit**: Fits a gaussian function onto the histogram of an IRF measruement (imaging a sample of gold nanoparticles)
* **Hist_monoexp_MLE_foreground**: Sums histograms of all pixels in the foreground of a FLIM image and fits a single exponential function to the data using MLE error estimation
* **Hist_monoexp_MLE_foreground_AllFolder**: Same as **Hist_monoexp_MLE_foreground** but repeats for all the FLIM images in a folder and saves results to a csv file.
* **Hist_monoexp_MLE_pixelwise**: Fits a single exponential function to the histogram of every pixel in an image data using MLE error estimation.
* **Hist_monoexp_MLE_pixelwise_AllFolder**: Same as **Hist_monoexp_MLE_pixelwise** but repeats for all the FLIM images in a folder

### PhasorDistribution
Display and analysis of phasor distributions of STED-FLIM images

* **PhasorDistribution_STED_EllipseFit_2Species** Analysis of phasor distribution trajectories as a function of STED depletion power by fitting ellipses on phasor distributions. Saves metrics (shortest distance between ellipses, IOU, distance between ellipse centroids , ...) to a csv file and graphs properties (mean phasor centroid, ellipse dimensions) for each depletion power.
<p align="center">    
<img src="images/2Dye_Centroids.png" width="270"/>
</p>

 <p> 
   <center><em>Mean ellipse and centroid for each depletion power </em></center>
 </p> 

* **ReadCSVtoGraphs_EllipseFit**: Reads csv file produced by **PhasorDistribution_STED_EllipseFit_2Species** and plots the different metrics as a function of depletion power
* **PhasorDistribution_STED_3DGraph** Makes 3D graph of phasor distributions of 2 stainings as a function of STED depletion power 
<p align="center">
<img src="images/3D_Phasor.png" width="400"/>
 <p> 
   <center><em>3D graph of phasors of 2 dyes as a function of depletion power</em></center>
 </p> 
</p>

* **PhasorDistribution_PrevsPostCalibration**: Plots a phasor before and after the calibration with the IRF
<p align="center">
<img src="images/PrevsPostCalibration.png" width="270"/>
 <p> 
   <center><em>phasors before and after application of calibration</em></center>
 </p> 
</p>

* **PhasorDistribution_RawvsMedianvsCWF**: Plots the same phasor with different filtering techniques (raw,median and CWF)
<p align="center">
<img src="images/RawvsFilters.png" width="400"/>
 <p> 
   <center><em>phasors before and after application of median and CWF filters</em></center>
 </p> 
</p>


### SPLIT-STED
1 Species [SPLIT-STED](https://doi.org/10.1039/C8NR07485B) to improve resolution of single staining STED-FLIM images
* **SPLIT-STED_CWF_multiSTEDpercent_Metrics**: 1 species SPLIT-STED of all images in a folder with [*CWF filtering*](https://doi.org/10.1364/BOE.420953) of phasors. Saves resolution, photobleaching and SQUIRREL metrics of the input and output images in a csv file
* **SPLIT-STED_Median_multiSTEDpercent_Metrics**: 1 species SPLIT-STED of all images in a folder with *median filtering* of phasors. Saves resolution, photobleaching and SQUIRREL metrics of the input and output images in a csv file
  
<img src="images/SPLIT_csv.png" width="640"/>  
<p float="middle">    
<img src="images/PhasorSPLIT_raw.png" width="320"/>
<img src="images/PhasorSPLIT.png" width="320"/>  
</p>
<p float="middle">    
<img src="images/Homer_STOrange_STEDPowerBleach_MediumAcq_1_10_20percentSTED_Confocal.png" width="210"/>
<img src="images/Homer_STOrange_STEDPowerBleach_MediumAcq_1_10_20percentSTED_STED.png" width="210"/>
<img src="images/Homer_STOrange_STEDPowerBleach_MediumAcq_1_10_20percentSTED_SPLIT_STED_Median_0.2.png" width="210"/>  
</p>
<img src="images/SPLIT_Fractmap.png" width="320"/>

* **ReadCSVstoGraphs_EllipsesvsSPLITresolution**: Read csv files generated by SPLIT-STED codes and **PhasorDistribution_STED_EllipseFit_2Species** and produces graphs


### Unmixing
Algorithms to separate dyes of different lifetimes in Confocal-FLIM and STED-FLIM images using [linear systems of equations in phasor space](https://doi.org/10.1088/2050-6120/ab8570)
* **Unmixing_2SpeciesSPLITSTED.ipynb** is a Jupyter Notebook implementation of 2 Species SPLIT-STED. A small set of example data is downloaded by the script from  [the paper website](https://s3.valeria.science/flclab-2-species-split-sted/index.html)into an **Example_data** subfolder 
* **Unmixing_2SpeciesConfocalFLIM__allfolder** : Separate phasors of Confocal-FLIM images of double-stained samples into two fractional components. Uses the same pair of control images for the entire folder of mixed images
* **Unmixing_2SpeciesSTEDFLIM_allfolder**: Separate phasors of STED-FLIM images of double-stained samples into two fractional components. Uses the same control images for all mixture images acquired with the same depletion power (1 pair of controls per depletion power)
* **Unmixing_2SpeciesSPLITSTED_allfolder** :  Separate phasors of STED-FLIM images of double-stained samples into three fractional components. Uses the same control images (4 per dye, confocal and 3 depletion powers) for all mixture images.

### Simulation
Creates simulated double-staining images by summing single-staining STED-FLIM images and performs different unmixing algorithms to evaluate their performance in a controlled setting.
* **HistogramMLE_MonoExponential_Simulationimages**: Generate simulated images and perform MLE histogram fitting to generate color-coded lifetime images of the ground truth and mixed images

<p float="middle">    
<img src="images/MLELifetime_IntensityComposite_0_Bassoon_CF594.png" width=210/>
<img src="images/MLELifetime_IntensityComposite_0_PSD95_STORANGE.png" width=210/>
<img src="images/MLELifetime_IntensityComposite_0_Mixture.png" width=210/>
</p>

* **Simulation_2SpeciesSTED-FLIM_looppowers** : Generate simulated images and perform 2 Species STED-FLIM unmixing. Calculate metrics (resolution, squirrel) and save them to a csv file. Also save input images and unmixing results (color-coded phasor and images) for each simulated combination.
* **Simulation_2SpeciesSPLIT-STED_looppowers**: Generate simulated images and perform 2 Species SPLIT-STED unmixing. Calculate metrics (resolution, squirrel) and save input images and unmixing results (color-coded phasor and images)
* **ReadCSVstoGraphs_Simulation** : Read the csv files produced by **Simulation_2SpeciesSTED-FLIM_looppowers** and **Simulation_2SpeciesSPLIT-STED_looppowers** and produce comparison graphs of the performance metrics of both unmixing methods as a function of the depletion power.

