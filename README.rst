Signature_PyProc
################

Version 0.0.1: *Under development.* 

Overview
--------

Post-processing and basic analysis of ice and ocean data from Nortek Signature
ADCPs. 

Designed for applications where the instrument is deployed looking upward below
the ocean surface.

  **GOAL:**
   
  **Easy and explicit post-processing to obtain scientific-quality estimates of**
    - *Sea ice draft* 
    - *Sea ice drift velocity*
    - *Upper ocean current velocity*
 

Inputs to the process are: 

- *.mat* files exported by Norteks `SignatureDeployment
  <https://www.nortekgroup.com/software>`_ software *(required)* 
- Time series of atmospheric pressure during the deployment *(strongly
  recommended)* 
- Time series of ocean temperature/salinity in the upper ocean *(recommended)* 
- Magnetic declination correction (single value or time series) *(necessary for
  correct velocity directions)* 

Core functionality
'''''''''''''''''''

**Read one or multiple (** *.mat* **) data files and transform to Python
dictionary**

- Reads and concatenates data from *.mat*-files produced using Nortek's
  `SignatureDeployment <https://www.nortekgroup.com/software>`_ software from
  *.ad2cp* data files output by the instrument. 
  
- Stores the data in an `xarray Dataset
  <https://docs.xarray.dev/en/stable/generated/xarray.Dataset.html>`_.

- Present version does not read *Burst* or *Waves* data - only *Average*
  (altimeter, ocean velocities) and *AverageIce* (ice velocities).

- The single *Average_Time* dimension is reshaped to two dimensions (*TIME*:
  time stamp of each ensemble, and *SAMPLE*: number of sample in ensemble). 

  - This is useful because the instrument usually samples in ensembles, e.g.
    collects 50 samples at 1 Hz once every 20 minutes. We tyopically want to do
    statistics on each ensembles to arrive at one value per *TIME*.      


**Estimate sea ice presence based on the** *Figure-of-Merit* **metric**

- Classifies each sample as ice/no ice based on the Figure-of-Merit (FOM) value
  of the ice velocity measurement of the four slanted ADCP beams. FOM is an
  indication of the Doppler noise of each ping; low values indicate ice and high
  values indicate open water.

- Sea ice presence is confirmed if FOM of all four beams is below a set
  threshold (default = 10 000). An estimate of "sea ice concentration" is
  calculated as the fraction of samples within an ensemble classified as ice.
  
  (*NOTE: This "sea ice concentration" is most meaningful when averaged over a
  longer time period, e.g. daily.*)

- Alternative estimates of sea ice presence/concentration (suffix *_FOM*) are
  made by a less conservative criterion: requiring that FOM of *at least one* of
  four beams is below the threshold. These are not recommended as they tend to
  give false positives for ice.


**Append data from external sources**

- Take for example a record of ocean temperature from another instrument, or sea
  ice concentration from remote sensing product, and add it to the sig500
  dataset interpolated onto the native time grid. This is useful for adding CTD
  variables (for sound speed corrections) or atmospheric pressure (for
  instrument depth correction), but can also be useful for analysis of sig500
  data in combination with e.g. remote sensing products. 

**Calculate ice draft based on altimeter data**

- Atmospheric pressure (from e.g. atmospheric reanalysis) and preferably ocean
  salinity and temperature (from e.g. a moored CTD sensor) should be supplied in
  order to get good quality results.
- Ice presence is determined based on the FOM (Figure-of-Merit) reading of the
  slanted beams. 
- *(Optional but recommended):* Perform an ad-hoc correction to the sound speed
  such that the long-term open water mean value is at 0 m.

**Basic editing/post-processing of ice and ocean velocity**

- Applies magnetic declination correction (rotating the velocity vector).
- Thresholds of FOM, quality, and speed used to filter out bad ice drift data.
- Various thresholds (speed, quality) used to filter out ocean velocity data. -
  Near- and above-surface measurements 

**Ensemble processing**

- Mean or median averaging to arrive at a single measurement per ensemble (per
  bin).
- Individual ensembles are retained only if they pass various quality tests and
  tests of internal consistency.

**Depth interpolation** - Interpolation fixed-depth.

**Analysis and visualization** - Pr

**Export** - Functionality to export to smaller "analysis" xarray Dataset where
unused
  variables are removed.
- Functionality to export to netCDF file (.nc). 

  *NOTE:* The resulting file *should* be formatted according to CF conventions.
  Do check your dataset closely before publication, however.

****



Modules
'''''''''''''''''''

´´sig_funcs.py´´

Functions for loading one or several *.mat* files from a deployment. Reads to
desired format, reshapes to the desired (TIME, SAMPLE) 2d shape, adds some
metadata, stores as an *xarray* Dataset.

Function for calculating tilt from pitch and roll.

Function for estimating ice presence/concentration. 

´´sig_append_external.py´´

Functions to append external datasets to an xarray Dataset containing Nortek
Signature data. 

- General function for adding and interpolating any time series data:

Some specialized wrapper functions used for loading data that needs to be
formatted correctly in later operations:

- Add CTD data and compute sound speed (for ice draft calculations)
- Add air pressure (for instrument depth corrections)
- Add magnetic declination (for correction of velocity directions)


Dependencies
-------------

*Signature_PyProc* is a Python package, and requires Python 3 (will not work on
1. and has currently only been tested on 3.8).

**Standard libraries:**

- ``numpy`` 
- ``scipy`` 
- ``matplotlib`` 

**Other:**

- `xarray <https://docs.xarray.dev/en/stable/>`_ - data are stored and
  manipulated as xarray *Dataset* objects.
- `GSW-Python <https://teos-10.github.io/GSW-Python/>`_ - used for computation
  of depth from pressure as well as density/sound speed/etc from CTD
  measurements.
 


Basic example
-------------

TBW

Version history
'''''''''''''''''''

Currently under development.