'''
Functions to append external datasets to an xarray Dataset containing Nortek
Signature data. 

- General function for adding and interpolating any time series data:

Some specialized wrapper functions used for loading data that 
needs to be formatted correctly in later operations:

- Add CTD data and compute sound speed (for ice draft calculations)
- Add air pressure (for instrument depth corrections)
- Add magnetic declination (for correction of velocity directions)

TO DO:
- Check air pressure wrapper.
- Think about whether the ERA-5 picker should be retained here or elsewhere
  (otherwise remove)

'''

import numpy as np
import xarray as xr
from scipy.interpolate import interp1d
import gsw
from sigpyproc.sig_funcs import mat_to_py_time
from matplotlib.dates import date2num

def add_to_sigdata(dx, data, time, name, attrs = None, time_mat = False):
    '''
    Adds a time series to the Signature dataset. Interpolates onto the "time"
    coordinate (one entry per ensemble).

    Used in the functions append_ctd() and append_slp(), but can also be 
    useful for appending e.g. remote sensing sea ice data for 
    comparison/validation.

    Inputs
    ------
    dx: xarray dataset with signature data. 
    data: Time series data
    time: Time grid of data (python epoch unless time_mat = True) 
    name: Name of the new variable (string) 
    attrs: Attributes of the new variable (dictionary). Good place to
           include "units", "long_name", etc..
    time_mat: Set to True if *time* is matlab epoch 
              (False/default: python default epoch)

    Outputs
    -------
    dx: The xarray dataset including the new variable.

    '''

    if time_mat:
        time = mat_to_py_time(time)

    # Interpolatant of the time series
    data_ip = interp1d(time, data, bounds_error = False)

    # Add interpolated data to dx
    dx[name] = (('TIME'), data_ip(dx.TIME.data), attrs)

    return dx


def append_ctd(dx, temp, sal, pres, CTDtime, instr_SN = None, instr_desc = None, 
                time_mat = False, pressure_offset = 0):
    '''
    Read data from a moored CTD - used for sound speed corrections etc. 
    Converts to TEOS-10 and computes sound speed using the gsw module. 

    Interpolates onto the *time* grid of the sig500 data. 

    Note: *temp, sal, pres* should be on the same *time* grid, 

    Inputs
    ------

    dx: xarray dataset with signature data.
    
    temp: In-situ temperature [C].
    salt: Practical salinity [].
    pres: Ocean pressure [dbar]. 
    CTDtime: Time stamp of CTD measurements .

    Outputs
    -------
    dx: The xarray dataset including the new SA, CT, pres_CTD, and 
        sound_speed_CTD variables.
    '''

    SA = gsw.SA_from_SP(sal, pres, dx.lon, dx.lat)
    CT = gsw.CT_from_t(SA, temp, pres)
    ss = gsw.sound_speed(SA, CT, pres)

    attrs_all = {'Instrument description': instr_desc, 'Instrument SN':None,
            'note':'Calculated using the gsw module. Linearly'
            ' interpolated onto Sig500 time grid.'}
    
    add_to_sigdata(dx, SA, CTDtime, 'SA_CTD', 
                attrs = {'long_name':'Absolute Salinity', 'units':'g kg-1',
                  **attrs_all},
                time_mat = time_mat)
    add_to_sigdata(dx, CT, CTDtime, 'CT_CTD', 
                attrs = {'long_name':'Conservative Temperature', 'units':'degC',
                    **attrs_all},
                time_mat = time_mat)
    add_to_sigdata(dx, CT, CTDtime, 'pres_CTD', 
                attrs = {'long_name':'Pressure (CTD measurements)', 
                    'units':'dbar', **attrs_all},
                time_mat = time_mat)
    add_to_sigdata(dx, ss, CTDtime, 'sound_speed_CTD', 
                attrs = {'long_name':'Sound speed', 'units':'m s-1',
                    **attrs_all},
                time_mat = time_mat)


def append_atm_pres(dx, slp, slptime, attrs = None, 
                    time_mat = False):
    '''
    Append sea level pressure from e.g. ERA-5. Note that the
    pressure units should be dbar.

    Interpolates onto the *time* grid of the sig500 data
    and adds to the sig500 data as the variable *SLP*. 

    Inputs
    ------

    dx: xarray dataset with signature data.
    
    slp: Sea level atmospheric pressure [dbar].
    slptime: Time stamp of slp.
    attrs: Attributes (dictionary).
    time_mat: Set to True if *slptime* is matlab epoch 
              (False/default: python default epoch).

    Outputs
    -------
    dx: The xarray dataset including the SLP variable.
    '''

    # Modify the attribute dictionary (specified "units" and "long_name"
    # overrides default ones).
    attrs_all = {'long_name':'Sea level pressure', 'units':'db'} 
    if attrs:
        for nm in attrs:
            attrs_all[nm] = attrs[nm]

    # Append to sig500 data
    add_to_sigdata(dx, slp, slptime, 'SLP', 
                attrs = attrs_all,
                time_mat = time_mat)


def append_magdec(dx, magdec, magdectime = False, attrs = None, 
                    time_mat = False):
    '''
    Append magnetic declination angle, used for correcting the heading of 
    observed velocities. 

    Magnetic declination can be supplied as a fixed number or as a 
    time-varying quantity (useful for longer deployments or cases
    where the position is not fixed).

    Appended to the sig500 data as the variable *magdec* - 
    either as a single number or interpolated onto the *time* grid
    of the sig500 data. 

    Inputs
    ------

    dx: xarray dataset with signature data.
    
    magdec: Magnetic declination angle [degrees] - single number or array
            of time-varying declination.
    magdectime: Time stamp of magdec if it is time-varying.
    attrs: Attributes (dictionary).
    time_mat: Set to True if *slptime* is matlab epoch 
              (False/default: python default epoch).

    Outputs
    -------
    dx: The xarray dataset including the magdec variable.
    '''

    # Modify the attribute dictionary (specified "units" and "long_name"
    # overrides default ones).
    attrs_all = {'long_name':'Magnetic declination', 'units':'degrees'} 
    if attrs:
        for nm in attrs:
            attrs_all[nm] = attrs[nm]

    # Append to sig500 data
    if hasattr(magdec, '__iter__'): # If this is an array of several 
                                    # magdec values
        if not magdectime:  
            raise Exception('Looks like you supplied a time-varying'
            '*magdec* but not the required time stamps *magdectime*..')
        else:
            add_to_sigdata(dx, magdec, magdectime, 'magdec', 
                attrs = attrs_all,
                time_mat = time_mat)

    else:
        dx['magdec'] = ((), magdec, attrs_all)

    return dx


# ERA-5 retrieval deprecated for now (decided for now that this is best suited elsewhere as users may 
# do this differently. This way is a bit cumbersome anyways..)

if False:
    def get_era5(dx, temp_2m = False, wind_10m = False):
        '''
        Read ERA-5 variables in the nearest model grid cell:

        - Sea level pressure [db] (Used for depth correction)
        - 2-m air temperature [deg C] (Optional - toggle with *temp_2m=True*)
        - 10-m wind components [m s-1] (Optional - toggle with *wind_10m=True*)

        Adds to the sig500 dictionary as 1D variables interpolated to the *time*
        grid.
        
        Accessing hourly ERA-5 data over OpenDap from the Asia-Pacific Data 
        Research Center (http://apdrc.soest.hawaii.edu/datadoc/ecmwf_ERA5.php).
        
        Note: Accessing the data can take as much as tens of minutes.
        '''
        # Loading remote datasets
        # This operation can take 10-20 seconds
        era5_url = ('http://apdrc.soest.hawaii.edu:80/dods/public_data/'
                    'Reanalysis_Data/ERA5/hourly/')

        print('Connecting to ERA5 MSL remote dataset..')
        era5_msl = xr.open_dataset(era5_url + 'Surface_pressure')
        if temp_2m:
            print('Connecting to ERA5 2-m air temperature remote dataset..')
            era5_msl = xr.open_dataset(era5_url + '2m_temperature')
        if wind_10m:
            print('Connecting to ERA5 10-m wind remote datasets..')
            era5_uwind = xr.open_dataset(era5_url + 'U_wind_component_10m')
            era5_vwind = xr.open_dataset(era5_url + 'V_wind_component_10m')
        print('..done.')

        # Find overlapping time stamps 
        t_era_full = date2num(era5.time.data)
        # Start and end points of sig500 deployment
        tind0_era = np.searchsorted(t_era_full, dx.time[0])-2
        tind1_era = np.searchsorted(t_era_full, dx.time[-1])+2
        tsl_era = slice(tind0_era, tind1_era)