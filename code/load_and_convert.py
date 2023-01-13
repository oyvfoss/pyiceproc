# Extract data from .mat files to useful format for further processing
# Modified on 24.08.22  to work when we don't have AverageIce/RawAltimeter data

from scipy.io import loadmat
from itertools import chain
#from oyv.data.signature import sig_funcs
import sig_funcs
import xarray as xr
import glob 

ddir_in = '/media/oyvindl/ratatosk/data/mooring_data/signature500/'
#ddir_out = '/media/sf_work/projects/sig500/data/proc/v1/'
ddir_out = '/media/oyvindl/ratatosk/data/mooring_data/signature500/'

mcoords = pickle_load('/home/oyvindl/work/arven_etter_nansen/data/aen_coords_v2.p')

#mf = loadmat(ddir + fn)
#bb = Bunch(b)
moorings = ['m1_1',  'm1_2', 'm2_1', 'm2_2', 'm2_3']
moorings = ['m1_3']

include_raw_altimeter = False

for mooring in moorings:

    outfn = '%s%s/proc/%s_signature_collected_v1.nc'%(ddir_out, mooring, mooring)
    
    # Get a list of mat files (not including avgd files)


    flist_all = glob.glob(ddir_in + '%s/converted/*.mat'%mooring)
    flist_a = np.sort([fn for fn in flist_all if 'avg' not in fn
                         and '_01_' not in fn and '_02_' not in fn])
    flist_b = np.sort([fn for fn in flist_all if 'avg' not in fn
                         and '_01_' in fn])
    flist_c = np.sort([fn for fn in flist_all if 'avg' not in fn
                         and '_02_' in fn])
    flist = np.concatenate((flist_b, flist_c, flist_a))

    pressure_offsets = np.array([])
   # data_arrays = []
    first = True
    for fn in flist:
        b = sig_funcs.sig_mat_to_dict(fn)


        # OBTAIN COORDINATES
        coords = {
            'time_average':time.mat_to_python(b['Average_Time']),
            'bins': np.arange(b['Average_VelEast'].shape[1]),
            'xyz': np.arange(3),
             }
        try: # If we have AverageRawAltimeter: Add this as well  
            coords += {'beams':np.arange(
                len(b['AverageRawAltimeter_BeamToChannelMapping'])),
                'along_altimeter':np.arange(
                    b['AverageRawAltimeter_AmpBeam5'].shape[1])}
        except:
            print('No *AverageRawAltimeter*')

        # CREATE AN XARRAY DATASET TO BE FILLED     
        dx = xr.Dataset(coords =coords)

        # Check whether we have AverageIce fields..
        try:
            dx['time_average_ice'] = (('time_average'), 
                time.mat_to_python(b['AverageIce_Time']))
        except:
            print('Did not find AverageIce..')

        # IMPORT DATA AND ADD TO XARRAY DATASET
        for key in b.keys():
            print('%s: %s..\r'%(fn[-10:], key), end = '') 
            

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

        pressure_offsets = np.append(pressure_offsets, b['conf']['PressureOffset'])


        if first: 
            DX = dx
            dx0 = dx.copy()
            first = False
        else:
            print('%s: CONCATENATING (print: first = %s)..\r'%(fn[-10:], first), end = '') 
            try:
                DX = xr.concat([DX, dx], dim = 'time_average')
            except:
                print('Failed at %s'%(fn[-10:]))


    DX.time_average.attrs['description'] = ('Time stamp for'
     ' "average" fields. Source field: *Average_Time*. Converted'
    ' using sig_functs.mat_to_python().') 
    if include_raw_altimeter:


        DX.time_raw_altimeter.attrs['description'] = ('Time stamp for'
            ' "AverageRawAltimeter" fields. Source field:'
            ' *AverageRawAltimeter_Time*. Converted'
            ' using time.mat_to_python.') 
 
    DX.bins.attrs['description'] = ('Number of velocity bins.') 
    DX.beams.attrs['description'] = ('Beam number (not 5th).') 
    DX.xyz.attrs['description'] = ('Spatial dimension.') 
    try: #(If we have this)
        DX.along_altimeter.attrs['description'] = ('Index along altimeter.') 
    except: pass

    units_str = ''
    for unit_ in b['units']:
        units_str += '%s: %s\n'%(unit_, b['units'][unit_])
    DX.attrs['units'] = units_str

    desc_str = ''
    for desc_ in b['desc']:
        desc_str += '%s: %s\n'%(desc_, b['desc'][desc_])
    DX.attrs['desc'] = desc_str

    conf_str = ''
    for conf_ in b['conf']:
       conf_str += '%s: %s\n'%(conf_, b['conf'][conf_])
    DX.attrs['conf'] = conf_str

    if len(np.unique(pressure_offsets)) == 1:
        DX.attrs['pressure_offset'] = pressure_offsets[0]
    else:
        DX.attrs['pressure_offset'] = pressure_offsets

    DX.attrs['lat'] = mcoords[mooring].lat
    DX.attrs['lon'] = mcoords[mooring].lon

    print('\nSAVING TO:  .. %s ..'%outfn[-40:])
    DX.to_netcdf(outfn, mode = 'w')
    print('..done. \n ')
    DX.close()
