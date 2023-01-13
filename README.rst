Signature_PyProc
################

Version 0.0.1: *Under development.* 

Overview
--------

Post-processing and basic analysis of ice and ocean data from Nortek Signature ADCPs.

Works with *.mat* files exported by Norteks *SignatureDeployment* software

Core functionality
'''''''''''''''''''

**Read one or multiple (** *.mat* **) data files and transform to Python dictionary.**

- The single time dimension is reshaped to two dimensions (*ensemble_time*,
   *ensemble_sample*). This is useful because this instruments usually samples in
   ensembles, e.g. collects 50 samples at 1 Hz once every 20 minutes. We often
   want to do statistics on each ensembles to arrive at one value per
   *ensemble_time*.      

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
  slanted beams. Further refining may be applied here eventually.. 
- Optional: Perform an ad-hoc correction to the sound speed such that the smoothed
  value in open water is at 0 m.

**Basic editing/post-processing of ice and ocean velocity**

- Applies magnetic declination correction (rotating the velocity vector).
- Thresholds of FOM and speed used to filter out bad ice drift data.
- Various thresholds used to filter out ocean velocity data.  

**Ensemble processing**

- Mean or median averaging to arrive at a single measurement per ensemble (per bin).
- Individual ensembles are retained only if they pass various quality tests and
  tests of internal consistency.

****



Modules
'''''''''''''''''''

´´sig_funcs.py´´

Functions for loading one or several *.mat* files from a deployment. Reads to
desired format, reshapes to the desired (TIME, SAMPLE) 2d shape, adds some metadata,
stores as an *xarray* Dataset.

Function for calculating tilt from pitch and roll.

´´sig_append_external.py´´

Functions to append external datasets to an xarray Dataset containing Nortek
Signature data. 

- General function for adding and interpolating any time series data:

Some specialized wrapper functions used for loading data that 
needs to be formatted correctly in later operations:

- Add CTD data and compute sound speed (for ice draft calculations)
- Add air pressure (for instrument depth corrections)
- Add magnetic declination (for correction of velocity directions)


Dependencies
-------------

*Signature_PyProc* is a Python package, and requires Python 3 (will not work on 2. 
and has currently only been tested on 3.8).

**Dependence, standard libraries:**

- ``numpy`` 
- ``scipy`` 
- ``matplotlib`` 

**Other:**

- `xarray <https://docs.xarray.dev/en/stable/>`_ - data are stored and manipulated 
  as xarray *Dataset* objects.
- `GSW-Python <https://teos-10.github.io/GSW-Python/>`_ - used for computation of
  sound speed.
 




**EVERYTHING BELOW IS COPY-PASTED FROM ELSEWHERE - IGNORE**





*Adcpyproc* is a Python package, and requires Python 3 (will not work on 2. 
and has currently only been tested on 3.8).

Beyond standard Python libraries such as *numpy*, *scipy*, *matplotlib*, etc, 
*adcpyproc* depends on 
- `pickle <https://docs.python.org/3/library/pickle.html>`_ (for saving python
  dictionaries).
- `geomag <https://pypi.org/project/geomag/>`_ (for computing time-and
  location-dependent magnetic declination based on the World Magnetic Model for
  compass corrections).  
- `netCDF4 <https://unidata.github.io/netcdf4-python/>`_ for exporting to
   netCDF.
 
To process raw *.000* files output from RDI ADCPs, you need to use RDI's WinADCP
software.

Installing
----------

To install, you should be able to do something like this:

1. Obtain a copy of *adcpyproc*. Either: 
    - Clone the git repository into a suitable location.
    - Download the zip file and unpack into a suitable folder.

2. Navigate to the top folder of the pproject (the one containing ``setup.py`` ).
3. If you haven't already: Install *pip*. For conda systems: ``conda install pip``.
4. Install *adcpyproc* using ``pip install -e .``.
    - The ``-e`` flag should be used while the code is still in the initial
      development phase.
5. You should now be able to load the module using, e.g.:
    - ``import adcpyproc``
    - ``from adcpyproc import rdi_adcp``

Example usage
-------------

Minimal example
+++++++++++++++

An example of processing steps which may be considered a minimum for ADCP
processing.


::

    # Import adcpyproc
    from adcpyproc import rdi_adcp

    # Set the path to the .mat file output from WinADCP
    fn_mat = 'path/to/matfile/file_from_winadcp.mat'

    # Load data from a WinADCP-produced .mat-file into an RdiObj object
    d = rdi_adcp.RdiObj(fn_mat)           

    # Chop away deployment and recovery
    d.remove_ship_time()   

    # Reject bins where the median depth is above the surface.
    d.reject_surface_bins()     

    # Set lat and lon and correct current direction for magnetic declination.
    d.set_latlon(80, 30)                   
    d.correct_magdec()   

    # MASKING/FLAGGING BASED ON THRESHOLDS 
    d.mask_umax(thr_uvamp=100)        # Mask entries where u or v exceed 100 cm/s amplitude.
    d.mask_pg(thr_pg=75)              # Mask entries where percent good (PG1 + PG4) is below 75%.
    d.mask_amp_a(thr_amp=64)          # Mask entries where beam mean amplitude is less than 64 db.
    d.mask_errvel(thr_errvel=50)      # Mask entries where error velocity is greater than 50 cm/s.
    d.mask_ca(thr_cor=30)             # Mask entries where the mean beam correlation in two or  
                                      # more beams is below 45 counts.
    d.mask_tilt(thr_tilt=20)          # Mask entries where tile exceeds 20 degrees.

    # Print a string showing a summary of the dataset (time/depth ranges, mean velocities, etc).
    d.print_summary()                      

    # Print a processing history listing the individual steps applied to the dataset.   
    d.print_proc()                                              

    # Save as matfile with the typically most important parameters (t, depth, u, v, ..)
    d.to_matfile('test_fn.mat', sparse = True)       


Extended example
+++++++++++++++++


::

    # Import adcpyproc
    from adcpyproc import rdi_adcp

    # Set the path to the .mat file output from WinADCP
    fn_mat = 'path/to/matfile/file_from_winadcp.mat'

    # Load data from a WinADCP-produced .mat-file into an RdiObj object
    d = rdi_adcp.RdiObj(fn_mat)           

    # -------------------------------------------------------------------------------------

    # Print some system parameters (instrument configuration)                 
    d.print_system_info()                  

    # Chop away deployment and recovery
    d.remove_ship_time()    

    # Adjust transducer and bin depths 3.2 m *downwards*               
    d.apply_depth_offset(3.2)              

    # Reject bins where the median depth is above the surface.
    d.reject_surface_bins()                

    # Reject the two bins nearest to transducer(rows 0 and 1)). Will prompt y/n.
    d.reject_rows([0, 1])                     

    # Set lat and lon and correct current direction for magnetic declination.
    d.set_latlon(80, 30)                   
    d.correct_magdec()                     

    # -------------------------------------------------------------------------------------

    ## Masking based on criteria (apply the relevant ones and modify the criteria) ##
    ## (Masks will end up as NaNs when exporting to matlab)

    d.mask_umax(thr_uvamp=100)           # Mask entries where u or v exceed 100 cm/s amplitude.
    # MASK_SURF_SIDELOBE NOT APPLIED YET!
    #d.mask_surf_sidelobe() (X)          # Mask entries falling within the estimated range of 
                                         # sidelobe interference of the surface.
    d.mask_pg(thr_pg=75)                 # Mask entries where percent good (PG1 + PG4) is below 75%.
    d.mask_amp_a(thr_amp=64)             # Mask entries where beam mean amplitude is less than 64 db.
    d.mask_errvel(thr_errvel=50)         # Mask entries where error velocity is greater than 50 cm/s.
    d.mask_ca(thr_cor=30)                # Mask entries where the mean beam correlation in two or  
                                         # more beams is below 45 counts.
    d.mask_ca_mean(thr_cor=30)           # Mask entries where the mean beam correlation is below
                                         # 45 counts.
    d.mask_w(thr_w=30)                   # Mask entries where the mean vertical is below 30 cm/s.
    d.mask_amp_jump(max_amp_increase=30) # Masking entries where the beam amplitude of any beam has 
                                         # a jump of 30 db or more (masking *after* the jump).
    d.mask_amp_jump(max_amp_increase=30, # Same, but also masks all entries *above* such jumps.
                    mask_above=True)
    d.mask_tilt(thr_tilt=20)             # Mask entries where tile exceeds 20 degrees.
    
    # -------------------------------------------------------------------------------------
    
    # calculate vertical shear (s2, shu, shv)
    d.calculate_shear()                   

    # Reject rows with less than 50% valid (unmasked) entries.
    d.reject_rows(masked_max=50)      

    # Print a string showing a summary of the dataset (time/depth ranges, mean velocities, etc).
    d.print_summary()                      

    # Print a processing history listing the individual steps applied to the dataset.   
    d.print_proc()                          

    # Export the dataset to a python Bunch.
    b = d.to_dict()                       

    # Save as matfile with the typically most important parameters (t, depth, u, v, ..)
    d.to_matfile('test_fn.mat', sparse=True)       

    # Save as pickled python dictionary (all parameters).
    d.to_matfile('test_fn_full.mat')      

     # Save as pickled python dictionary (all parameters).
    d.to_pickle('test_fn.p')              
    
    # Save as netcdf file (all parameters). TO_NETCDF4 NOT APPLIED YET!
    #d.to_netcdf('test_fn.nc')            