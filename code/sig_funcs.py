'''
Various functions for loading and concatenating Nortek Signature matfiles
produced by Nortek SignatureDeployment.

TO DO:

- Checking, fixing
- Does the *sig_mat_to_dict_join* do anything now?
- Make it easier/clearer to preserve/skip IBurst, RawAltimeter etc..
- Time: Specify epoch?

'''

##############################################################################

# IMPORTS

import numpy as np
from scipy.io import loadmat
import xarray as xr
import glob 
#
import pdb
##############################################################################

def matfiles_to_dataset(file_list, reshape = True, lat = None, lon = None,
                include_raw_altimeter = False, ):
    '''
    Read, convert, and concatenate .mat files exported from SignatureDeployment.

    Inputs:
    -------

    file_list: list of .mat files.
    fold: If True, reshape all time series from a single 'time_average'
          dimension to 2D ('TIME', 'SAMPLE') where we TIME is the 
          mean time of each ensemble and SAMPLE is each sample in the 
          ensemble. 
    lat, lon: Lat/lon of deployment (single point)
    include_raw_altimeter: Include raw altimeter signal if available.
                           (Typically on a single time grid) 

    Output:
    -------

    xarray Dataset containing the data.

    '''

    ###########################################################################
    # LOAD AND CONCATENATE DATA 

    first = True
    pressure_offsets = np.array([])

    if len(file_list)==0:
        raise Exception(
            'The *file_list* given to the function '
            'sig_funcs.matfiles_to_dataset() is empty.')

    for filename in file_list:
       
        dx, pressure_offset = _matfile_to_dataset(filename,
            lat = lat, lon = lon, 
            include_raw_altimeter = include_raw_altimeter)
        
        pressure_offsets = np.append(pressure_offsets, pressure_offset)

        if first: 
            DX = dx
            dx0 = dx.copy()
            first = False
        else:
            print('%s: CONCATENATING (print: first = %s)..\r'%(
                filename[-10:], first), end = '') 
            try:
                DX = xr.concat([DX, dx], dim = 'time_average')
            except:
                print('Failed at %s'%(fn[-10:]))


    ###########################################################################

    if len(np.unique(pressure_offsets)) == 1:
        DX.attrs['pressure_offset'] = pressure_offsets[0]
    else:
        DX.attrs['pressure_offset'] = pressure_offsets

    DX.attrs['lat'] = lat
    DX.attrs['lon'] = lon


    # Add tilt (from pitch/roll)
    DX = _add_tilt(DX)

    # Reshape
    if reshape:
        DX = reshape_ensembles(DX)

    return DX
##############################################################################

def reshape_ensembles(DX):
    '''
    Reshape all time series from a single 'time_average'
    dimension to 2D ('TIME', 'SAMPLE') where we TIME is the 
    mean time of each ensemble and SAMPLE is each sample in the 
    ensemble. 
    '''

    ###########################################################################
    # ADD A "time" COORDINATE (ONE ENTRY PER ENSEMBLE)

    Nt = len(DX.time_average)

    # Find indices where there is a time jump > 7 minutes
    time_jump_inds = np.where(np.diff(DX.time_average)*24*60>7)[0]

    # Start and end times of each ensemble
    ens_starts = np.concatenate([np.array([0]), time_jump_inds+1]) 
    ens_ends = np.concatenate([time_jump_inds, 
                           np.array([Nt-1])])

    # Use the mean time of each ensemble 
    t_ens = 0.5*(DX.time_average.data[ens_starts]
               + DX.time_average.data[ens_ends])
    Nens = len(t_ens)
    Nsamp_per_ens = int(Nt/Nens)
    print('%i time points, %i ensembles. Sample per ensemble: %i'%(
         Nt, Nens, Nsamp_per_ens))


    # NSW XARRAY DATASET  

    # Inheriting dimensions except time_average
    rsh_coords = dict(DX.coords)
    rsh_coords.pop('time_average') 

    # New coordinates: "TIME" (time stamp of each ensemble) and "SAMPLE" 
    # (number of samples within ensemble)

    rsh_coords['TIME'] = (['TIME'], t_ens, 
              {'units':'Days since 1970-01-01', 
               'long_name':('Time stamp of the ensemble averaged'
               ' measurement')})
    rsh_coords['SAMPLE'] = (['SAMPLE'], 
                np.int_(np.arange(1, Nsamp_per_ens+1)),
                {'units':'Sample number', 
               'long_name':('Sample number in ensemble '
               '(%i samples per ensemble)'%Nsamp_per_ens)})

    DXrsh = xr.Dataset(coords = rsh_coords)

    # Inherit attributes
    DXrsh.attrs = DX.attrs

    # Loop through variables, reshape where necessary
    for var_ in DX.variables:
        if DX[var_].dims == ('time_average',):
           # pdb.set_trace()
            DXrsh[var_] = (('TIME', 'SAMPLE'), 
                    np.ma.reshape(DX[var_], (Nens, Nsamp_per_ens)),
                    DX[var_].attrs)
        elif DX[var_].dims == ('bins', 'time_average'):

            DXrsh[var_] = (('bins', 'TIME', 'SAMPLE'), 
                    np.ma.reshape(DX[var_], (DX.dims['bins'], 
                        Nens, Nsamp_per_ens)),
                    DX[var_].attrs)
        elif DX[var_].dims == ('time_average', 'xyz'):
            DXrsh[var_] = (('TIME', 'SAMPLE', 'xyz'), 
                    np.ma.reshape(DX[var_], (
                        Nens, Nsamp_per_ens, DX.dims['xyz'])),
                    DX[var_].attrs)
        elif DX[var_].dims == ('time_average', 'beams'):
            DXrsh[var_] = (('TIME', 'SAMPLE', 'beams'), 
                    np.ma.reshape(DX[var_], (
                        Nens, Nsamp_per_ens, DX.dims['beams'])),
                    DX[var_].attrs)

    return DXrsh

##############################################################################

def _matfile_to_dataset(filename, lat = None, lon = None,
                include_raw_altimeter = False):
    '''
    Read and convert single .mat file exported from SignatureDeployment.

    Wrapped into *matfiles_to_dataset*.


    Inputs:
    -------

    file_list: list of .mat files.
    lat, lon: Lat/lon of deployment (single point)
    include_raw_altimeter: Include raw altimeter signal if available.
                           (Typically on a single time grid) 

    Output:
    -------
    dx: xarray Dataset containing the data.
    pressure_offset: Pressure offset used in the data.
    '''

    b = _sig_mat_to_dict(filename)

    # OBTAIN COORDINATES
    coords = {
        'time_average':mat_to_py_time(b['Average_Time']),
        'bins': np.arange(b['Average_VelEast'].shape[1]),
        'xyz': np.arange(3),
            }

    if include_raw_altimeter:
        try: # If we have AverageRawAltimeter: Add this as well  
            coords += {'beams':np.arange(
                len(b['AverageRawAltimeter_BeamToChannelMapping'])),
                'along_altimeter':np.arange(
                    b['AverageRawAltimeter_AmpBeam5'].shape[1])}
        except:
            print('No *AverageRawAltimeter*')

    # CREATE AN XARRAY DATASET TO BE FILLED     
    # Further down: concatenating into a joined dataset.
    dx = xr.Dataset(coords = coords)

    # Check whether we have AverageIce fields..
    try:
        dx['time_average_ice'] = (('time_average'), 
            mat_to_py_time(b['AverageIce_Time']))
    except:
        print('Did not find AverageIce..')

    # IMPORT DATA AND ADD TO XARRAY DATASET
    # Looping through variables. Assigning to dataset 
    # as fields with the appropriate dimensions.

    for key in b.keys():
        print('%s: %s..\r'%(filename[-10:], key), end = '') 
        
        # AverageIce fields
        if 'AverageIce_' in key:
            if b[key].ndim==0:
                dx[key] = ((), b[key])
            elif b[key].ndim==1:
                dx[key] = (('time_average'), b[key])
            elif b[key].ndim==2:
                dx[key] = (('time_average', 'xyz'), b[key])
    
        # Average fields
        elif 'Average_' in key:
            if b[key].ndim==0:
                dx[key] = ((), b[key])
            elif b[key].ndim==1:
                if len(b[key]) == dx.dims['time_average']:
                    dx[key] = (('time_average'), b[key])
                else:
                    dx[key] = (('beams'), b[key])
            elif b[key].ndim==2:
                if b[key].shape[1] == dx.dims['xyz']:
                    dx[key] = (('time_average', 'xyz'), b[key])
                elif b[key].shape[1] == dx.dims['bins']:
                    dx[key] = (('bins', 'time_average'), 
                            b[key].T)

        # AverageRawAltimeter fields
        elif 'AverageRawAltimeter_' in key:
            if include_raw_altimeter:
                if b[key].ndim==0:
                    dx[key] = ((), b[key])
                    
                elif b[dxkey].ndim==1:
                    if len(b[key]) == dx.dims['beams']:
                        dx[key] = (('beams'), b[key])
                    else:
                        dx[key] = (('time_raw_altimeter'), b[key])
                elif b[key].ndim==2:
                    if b[key].shape[1] == dx.dims['xyz']:
                        dx[key] = (('time_raw_altimeter', 'xyz'), b[key])
                    elif b[key].shape[1] == dx.dims['along_altimeter']:
                        dx[key] = (('along_altimeter', 
                                    'time_raw_altimeter'), b[key].T)


    # ASSIGN METADATA ATTRIBUTES
    # Units, description etc
    dx.time_average.attrs['description'] = ('Time stamp for'
     ' "average" fields. Source field: *Average_Time*. Converted'
    ' using sig_funcs.mat_to_py_time().') 
    dx.bins.attrs['description'] = ('Number of velocity bins.') 
    dx.beams.attrs['description'] = ('Beam number (not 5th).') 
    dx.xyz.attrs['description'] = ('Spatial dimension.') 
    
    if include_raw_altimeter:
        dx.time_raw_altimeter.attrs['description'] = ('Time stamp for'
            ' "AverageRawAltimeter" fields. Source field:'
            ' *AverageRawAltimeter_Time*. Converted'
            ' using mat_to_py_time.') 
        dx.along_altimeter.attrs['description'] = ('Index along altimeter.') 
 
    units_str = ''
    for unit_ in b['units']:
        units_str += '%s: %s\n'%(unit_, b['units'][unit_])
    dx.attrs['units'] = units_str

    desc_str = ''
    for desc_ in b['desc']:
        desc_str += '%s: %s\n'%(desc_, b['desc'][desc_])
    dx.attrs['desc'] = desc_str

    conf_str = ''
    for conf_ in b['conf']:
       conf_str += '%s: %s\n'%(conf_, b['conf'][conf_])
    dx.attrs['conf'] = conf_str

    # Read pressure offset

    pressure_offset = b['conf']['PressureOffset']

    return dx, pressure_offset


##############################################################################

def _sig_mat_to_dict(matfn, include_metadata = True, squeeze_identical = True,
                    skip_Burst = True, skip_IBurst = True, 
                    skip_Average = False,  skip_AverageRawAltimeter = False, 
                    skip_AverageIce = False, 
                    skip_fields = []):
    '''
    Reads matfile produced by SignatureDeployment to numpy dictionary.

    include_metadata = False will skip the config, description, and units
    fields. 
    
    squeeze_identical = True will make data arrays of identical entries into
    one single entry - e.g. we just need one entry, [N], showing the number of
    bins, not one entrye per time stamp [N, N, ..., N]. 
    '''

    # Load the mat file
    dmat = loadmat(matfn)

    # Create a dictionary that we will fill in below
    d = {}

    # Add the contents of the non-data fields
    for indict_key, outdict_key in zip(['Config', 'Descriptions', 'Units'],
                                ['conf', 'desc', 'units']):
        outdict = {}
        for varkey in dmat[indict_key].dtype.names:
            outdict[varkey] = _unpack_nested(dmat[indict_key][varkey])
        d[outdict_key] = outdict

    # Making a list of strings based on the *skip_* booleans in the function
    # call. *startsrtings_skip* contains a start strings. If a variable starts
    # with any of these srtings, we won't include it.
    
    startstrings = ['Burst_', 'IBurst_', 'Average_', 
                'AverageRawAltimeter_', 'AverageIce_', ]
    startstrings_skip_bool = [skip_Burst, skip_IBurst, skip_Average,
                skip_AverageRawAltimeter, skip_AverageIce,]
    startstrings_skip = tuple([str_ for (str_, bool_) in zip(startstrings, 
                        startstrings_skip_bool) if bool_])

    # Add the data. Masking NaNs and squeezing/unpacking unused dimensions
    for varkey in dmat['Data'].dtype.names:
        if (not varkey.startswith(startstrings_skip) 
           and varkey not in skip_fields):
           
            d[varkey] = np.ma.squeeze(np.ma.masked_invalid(
                    _unpack_nested(dmat['Data'][varkey])))
            
            if (d[varkey] == d[varkey][0]).all() and squeeze_identical:
                d[varkey] = d[varkey][0]

    # Return the dictionary
    return d

##############################################################################

if False: # DEPRECATED - NOW DOING THIS IN XARRAY?
    def sig_mat_to_dict_join(flist, **kwargs):
        '''
        E.g. skip_AverageIce = True
        '''
        D = sig_mat_to_dict(flist[0], include_metadata = True, 
                            **kwargs)
        sn0 = D['conf']['SerialNo']

        for nn, fn in enumerate(flist[1:]):
            print('Reading file %i/%i..\r'%(nn+1, len(flist)-1), end = '')
            d_ = sig_mat_to_dict(fn, include_metadata = True, 
                            **kwargs)

            if d['conf']['SerialNo'] != sn0:
                raise Exception('Instrument serial number changed. Should '
                    'only join files that are part of the same deployment!')

        d['conf']['SerialNo']

##############################################################################

def _unpack_nested(val):
    '''
    Unpack nested list/arrays down to the base level.
    Need this because data in the SignatureDeployment matfiles end up in
    weird nested structures.
    
    E.g. [['s']] -> 's'
    E.g. [[[[[1, 2], [3, 4]]]] -> [[1, 2], [3, 4]]
    '''
    
    unpack = True
    while unpack == True:
        if hasattr(val, '__iter__') and len(val) == 1:
            val = val[0]
            if isinstance(val, str):
                unpack = False
        else:
            unpack = False
    return val

##############################################################################


def mat_to_py_time(mattime):
    '''
    Convert matlab datenum (days) to Matplotlib dates (days).

    MATLAB base: 00-Jan-0000
    Matplotlib base: 01-Jan-1970
    '''

    mpltime = mattime - 719529.0

    return mpltime

##############################################################################

def _add_tilt(d):
    '''
    Calculate tilt from pitch/roll components. See Mantovanelli et al 2014 and
    Woodgate et al 2011.

    Input: xarray dataset with pitch and roll fields as read by
    matfiles_to_dataset
            
    '''

    tilt_attrs = {
        'units':'degrees', 
        'desc': ('Tilt calculated from pitch+roll'),
        'note':('Calculated using the function sig_funcs._add_tilt(). '
            'See Mantovanelli et al 2014 and Woodgate et al 2011.'),}

    try:
        d['tilt_Average'] = (('time_average'), 
            180 / np.pi* np.arccos(np.sqrt(
            1 - np.sin(d.Average_Pitch.data/180*np.pi) ** 2
            - np.sin(d.Average_Roll.data/180*np.pi) ** 2)))
        d['tilt_Average'].attrs  =  tilt_attrs
    except:
        print('(No AverageIce variables)')

    try:
        d['tilt_AverageIce'] = (('time_average'), 
            180 / np.pi* np.arccos(np.sqrt(
            1 - np.sin(d.AverageIce_Pitch.data/180*np.pi) ** 2
            - np.sin(d.AverageIce_Roll.data/180*np.pi) ** 2)))
        d['tilt_AverageIce'].attrs  = tilt_attrs
    except:
        print('(No AverageIce variables)')

    return d



