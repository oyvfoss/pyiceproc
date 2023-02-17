'''
SIG_VEL.PY

Functions for processing ocean and sea ice drift velocity.
'''

import numpy as np
import xarray as xr


def calculate_ice_vel(DX, avg_method = 'median'):
    '''
    Calculate sea ice drift velocity from the AverageIce 
    '''

    DX['uice'] = (('TIME', 'SAMPLE'), DX.AverageIce_VelEast.data, 
        {'units':'m s-1', 'long_name':'Eastward sea ice drift velocity',
        'details': 'All average mode samples'})
    DX['vice'] = (('TIME', 'SAMPLE'), DX.AverageIce_VelNorth.data, 
        {'units':'m s-1', 'long_name':'Northward sea ice drift velocity',
        'details': 'All average mode samples'})

    DX['uice'] = DX['uice'].where(DX.ICE_IN_SAMPLE)
    DX['vice'] = DX['vice'].where(DX.ICE_IN_SAMPLE)

    for key in ['uice', 'vice']:
        DX[key].attrs['processing_history'] = (
            'Loaded from AverageIce_VelEast/AverageIce_VelEast fields.\n')

    DX = DX._calculate_uvice_avg(DX, avg_method = avg_method)

    return DX 



def _calculate_uvice_avg(DX, avg_method = 'median'):
    '''
    Calculate ensemble average ocean velocity
    '''

    if avg_method=='median':
        DX['Uice'] = DX['uice'].median(dim = 'SAMPLE')
        DX['Vice'] = DX['vice'].median(dim = 'SAMPLE')
    elif avg_method=='mean':
        DX['Uice'] = DX['uice'].mean(dim = 'SAMPLE')
        DX['Vice'] = DX['vice'].mean(dim = 'SAMPLE')
    else:
        raise Exception('Invalid "avg_method" ("%s"). '%avg_method
         + 'Must be "mean" or "median".')

    DX.Uice.attrs = {'units':'m s-1', 
        'long_name':'Eastward sea ice drift velocity', 
        'details': 'Ensemble average (%s)'%avg_method, 
        'processing_history' : DX.uice.processing_history}  
    DX.Vice.attrs = {'units':'m s-1', 
        'long_name':'Northward sea ice drift velocity', 
        'details': 'Ensemble average (%s)'%avg_method, 
        'processing_history' : DX.vice.processing_history} 


    DX['Uice_SD'] = DX['uice'].std(dim = 'SAMPLE')
    DX['Vice_SD'] = DX['vice'].std(dim = 'SAMPLE')
    
    DX.Uice_SD.attrs = {'units':'m s-1', 
        'long_name':('Ensemble standard deviation of '
                'eastward sea ice drift velocity'), } 
    DX.Vice_SD.attrs = {'units':'m s-1', 
        'long_name':('Ensemble standard deviation of '
                'northward sea ice drift velocity'), }
    
    return DX

def calculate_ocean_vel(DX, avg_method = 'median'):
    '''
    Calculate sea ice drift velocity from the Average_VelEast/Average_VelNorth
    fields.
    '''
    # Calculate bin depths
    DX = _calculate_bin_depths(DX)

    # Extract u, v, data
    DX['uocean'] = (('BINS', 'TIME', 'SAMPLE'), DX.Average_VelEast.data, 
        {'units':'m s-1', 'long_name':'Eastward ocean velocity',
        'details': 'All average mode samples'})
    DX['vocean'] = (('BINS', 'TIME', 'SAMPLE'), DX.Average_VelNorth.data, 
        {'units':'m s-1', 'long_name':'Northward ocean velocity',
        'details': 'All average mode samples'})

    for key in ['uocean', 'vocean']:
        DX[key].attrs['details'] = 'All average mode samples'
        DX[key].attrs['units'] = 'm s-1'
        DX[key].attrs['processing_history'] = (
            'Loaded from Average_VelEast/Average_VelEast fields.\n')

    # Calculate sample averages
    DX = _calculate_uvocean_avg(DX, avg_method = avg_method )

    return DX 




def _calculate_uvocean_avg(DX, avg_method = 'median'):

    if avg_method=='median':
        DX['Uocean'] = DX.uocean.median(dim = 'SAMPLE')
        DX['Vocean'] = DX.vocean.median(dim = 'SAMPLE')
    elif avg_method=='mean':
        DX['Uocean'] = DX.uocean.mean(dim = 'SAMPLE')
        DX['Vocean'] = DX.vocean.mean(dim = 'SAMPLE')
    else:
        raise Exception('Invalid "avg_method" ("%s"). '%avg_method
         + 'Must be "mean" or "median".')

    DX.Uocean.attrs = {'units':'m s-1', 'long_name':'Eastward ocean velocity', 
        'details': 'Ensemble average (%s)'%avg_method, 
        'processing_history':DX.uocean.processing_history} 
    DX.Vocean.attrs = {'units':'m s-1', 'long_name':'Northward ocean velocity', 
        'details': 'Ensemble average (%s)'%avg_method, 
        'processing_history':DX.vocean.processing_history} 

    return DX


def _calculate_bin_depths(DX):
    '''
    Calculate time-varying depth of ocean velocity bins.

    From Nortek doc "Signature Principles of Operation":
        n-th cell is centered at a vertical distance from the transducer
        equal to: Center of n'th cell = Blanking + n*cell size

    '''

    dist_from_transducer = (DX.blanking_distance_oceanvel 
            + DX.cell_size_oceanvel*(1+np.arange(DX.N_cells_oceanvel)))

    DX['bin_depth'] = (DX.depth.mean(dim='SAMPLE').expand_dims(
                                dim = {'BINS':DX.dims['BINS']}) 
                       - dist_from_transducer[:, np.newaxis])
    DX['bin_depth'].attrs = {
            'long_name':'Sample-average depth of velocity bins',
            'units':'m', 
            'note':('Calculated as:\n\n'
            '   bin_depth = instr_depth - (blanking depth + n*cell size)\n\n'
            'where *n* is bin number and *inst_depth* the (sample-mean) depth'
            ' of the transducer.')}      

    return DX



def _uvoc_mask_sidelobe(DX):
    pass

def _uvoc_mask_range(DX, uv_max = 1.5, tilt_max = 5, 
        sspd_range = (1400, 1560), cor_min = 60,
        amp_range = (30, 85), max_amp_increase = 20):
    '''
    Set instances of 

    Setting uocean, vocean data points to NaN where: 

    1. Speed exceeds *uv_max* (m/s).
    2. Tilt exceeds *tilt_max* (degrees).
    3. Instrument-recorded sound speed is outside *sspd_range* (ms-1).
    4. Any one of the beam correlations is below *cor_min* (percent).
    5. Any one of the beam amplitudes are outside *amp_range* (db).
    6. There is a bin-to-bin amplitude jump of *max_amp_increase* (db)
       in any of the beams.
       - Also masking all beams *above* the jump. 

    '''
    # N_variables used for counting the effect of each step.
    N_start = np.float(np.sum(~np.isnan(DX.uocean)).data)

    # Create DX_uv; a copy of DX containing only uocean, vocean.
    # Then feeding these back into DX before returning.
    # (This is because we dont want the DX.where() operation
    # to affect other fields/expand dimensions unrelated to ocean 
    # velocities)
    DX_uv = DX[['uocean', 'vocean']]

    # Speed test
    DX_uv = DX_uv.where((DX.uocean**2 + DX.vocean**2)<uv_max**2)
    N_speed = np.float(np.sum(~np.isnan(DX_uv.uocean)).data)

    # Tilt test
    DX_uv = DX_uv.where(DX.tilt_Average<tilt_max)
    N_tilt = np.float(np.sum(~np.isnan(DX_uv.uocean)).data)
    
    # Sound speed test
    DX_uv = DX_uv.where((DX.Average_Soundspeed>sspd_range[0]) 
        & (DX.Average_Soundspeed<sspd_range[1]))
    N_sspd = np.float(np.sum(~np.isnan(DX_uv.uocean)).data)

    # Correlation test
    DX_uv = DX_uv.where((DX.Average_CorBeam1>cor_min) 
                | (DX.Average_CorBeam2>cor_min)
                | (DX.Average_CorBeam3>cor_min)
                | (DX.Average_CorBeam4>cor_min))
    N_cor = np.float(np.sum(~np.isnan(DX_uv.uocean)).data)

    # Amplitude test
    # Lower bound
    DX_uv = DX_uv.where((DX.Average_AmpBeam1>amp_range[0]) 
                | (DX.Average_AmpBeam2>amp_range[0])
                | (DX.Average_AmpBeam3>amp_range[0])
                | (DX.Average_AmpBeam4>amp_range[0]))
    # Upper bound
    DX_uv = DX_uv.where((DX.Average_AmpBeam1<amp_range[1]) 
                | (DX.Average_AmpBeam2<amp_range[1])
                | (DX.Average_AmpBeam3<amp_range[1])
                | (DX.Average_AmpBeam4<amp_range[1]))

    N_amp = np.float(np.sum(~np.isnan(DX_uv.uocean)).data)

    # Amplitude bump test

    # Find bumps from *diff* in the BIN S dimension
    is_bump = (
        (DX.Average_AmpBeam1.diff(dim = 'BINS')>max_amp_increase)
        | (DX.Average_AmpBeam2.diff(dim = 'BINS')>max_amp_increase)
        | (DX.Average_AmpBeam3.diff(dim = 'BINS')>max_amp_increase)
        | (DX.Average_AmpBeam4.diff(dim = 'BINS')>max_amp_increase))
    
    # Create a boolean (*True* above bumps)
    zeros_firstbin = xr.zeros_like(DX.uocean.isel(BINS=0))
    NOT_ABOVE_BUMP = xr.concat([zeros_firstbin, is_bump.cumsum(axis = 0)>0], 
                        dim = ('BINS'))<1
    DX_uv = DX_uv.where(NOT_ABOVE_BUMP)
    N_amp_bump = np.float(np.sum(~np.isnan(DX_uv.uocean)).data)


    proc_string = ('\nTHRESHOLD_BASED DATA CLEANING : '
     + '\nStart: %i initial valid samples.\n'%N_start
     + 'Dropping (NaNing samples where):\n'
     + '- # Speed < %.2f ms-1 # -> Dropped %i pts (%.2f%%)\n'%(
        uv_max, N_start-N_speed, (N_start-N_speed)/N_start*100)
     + '- # Tilt < %.2f deg # -> Dropped %i pts (%.2f%%)\n'%(
        tilt_max, N_speed-N_tilt, (N_speed-N_tilt)/N_speed*100)       
     + '- # Sound sp in [%.0f, %.0f] ms-1 # -> Dropped %i pts (%.2f%%)\n'%(
        *sspd_range, N_tilt-N_sspd, (N_tilt-N_sspd)/N_tilt*100)    
     + '- # Corr (all beams) < %.1f %% # -> Dropped %i pts (%.2f%%)\n'%(
        cor_min, N_sspd-N_cor, (N_sspd-N_cor)/N_sspd*100) 
     + '- # Amp (all beams) in [%.0f, %.0f] db # -> Dropped %i pts (%.2f%%)\n'%(
        *amp_range, N_cor-N_amp, (N_cor-N_amp)/N_cor*100)    
     + '- # Above amp bumps > %.0f db # -> Dropped %i pts (%.2f%%)\n'%(
        max_amp_increase, N_amp-N_amp_bump, (N_amp-N_amp_bump)/N_amp*100)   
     + 'End: %i valid samples.\n'%N_amp_bump 
        )

    for key in ['uocean', 'vocean']:
        DX[key] = DX_uv[key] 
        DX[key].attrs['processing_history'] += proc_string

    return DX



def rotate_vels_magdec(DX):
    '''
    Rotate ocean and ice velocities to account for magnetic declination.

    Only rotates processed fields (uocean, vocean, uice, vice) -
    not raw variables (Average_VelNorth, AverageIce_VelNorth, etc..).
    '''
    
    assert hasattr(DX, 'magdec'), ("Didn't find magnetic declination (no" 
        " *magdec* attribute). Run sig_append.append_magdec()..")

    DX0 = DX.copy()

    # Convert to radians
    magdec_rad = DX.magdec * np.pi/180

    # Make a documentation string (to add as attribute)
    magdec_mean = DX.magdec.mean().data
    magdec_str = ('Rotated CW by an average of %.2f degrees'%magdec_mean
                  +  ' to correct for magnetic declination. ')

    # Loop through different (u, v) variable pairs and rotate them
    uvpairs = [('uice', 'vice'), ('uocean', 'vocean'),] #
              #  ('Uice', 'Vice'), ('Uocean', 'Vocean')]

    uvstrs = ''

    for uvpair in uvpairs:

        if hasattr(DX, uvpair[0]) and hasattr(DX, uvpair[1]):
            uvc0_ = DX[uvpair[0]] + 1j* DX[uvpair[1]]
            uvc_rot = uvc0_ * np.exp(-1j*magdec_rad)

            DX[uvpair[0]].data = uvc_rot.real
            DX[uvpair[1]].data = uvc_rot.imag
            for key in uvpair:
                DX[key].attrs['processing_history'] += magdec_str+'\n'
            uvstrs += '\n - (%s, %s)'%uvpair

    if hasattr(DX, 'declination_correction'):
        inp_yn = float(input('Declination correction rotation has been '
        + 'applied to something before. \n -> Continue '
        + '(1) or skip new correction (0): '))

        if inp_yn ==1:
            print('-> Applying new correction.')
            DX.attrs['declination_correction'] = (
                '!! NOTE !! Magnetic declination correction has been applied more'
                ' than once - !! CAREFUL !!\n' 
                + DX.attrs['declination_correction'])
        else:
            print('-> NOT applying new correction.')
            return DX0
    else:
        DX.attrs['declination_correction'] = (
          'Magdec declination correction rotation applied to: %s'%uvstrs)
    try:
        DX = _calculate_uvocean_avg(DX, avg_method = 'median') 
    except:pass
    try:
        DX = _calculate_uvice_avg(DX, avg_method = 'median') 
    except:pass

    return DX

