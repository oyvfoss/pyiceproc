'''
SIG_ICE_VEL.PY

Functions for processing sea ice drift velocity.
'''
import numpy as np

def calculate_drift(DX, avg_method = 'median'):
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

    if avg_method=='median':
        DX['Uice'] = DX['uice'].median(dim = 'SAMPLE')
        DX['Vice'] = DX['vice'].median(dim = 'SAMPLE')
    elif avg_method=='mean':
        DX['Uice'] = DX['uice'].mean(dim = 'SAMPLE')
        DX['Vice'] = DX['vice'].mean(dim = 'SAMPLE')
    else:
        raise Exception('Invalid "avg_method" ("%s"). '%avg_method
         + 'Must be "mean" or "median".')

    DX.Uice.attrs = {'units':'m s-1', 'long_name':'Eastward sea ice drift velocity', 
        'details': 'Ensemble average (%s)'%avg_method} 
    DX.Vice.attrs = {'units':'m s-1', 'long_name':'Northward sea ice drift velocity', 
        'details': 'Ensemble average (%s)'%avg_method} 


    DX['Uice_SD'] = DX['uice'].std(dim = 'SAMPLE')
    DX['Vice_SD'] = DX['vice'].std(dim = 'SAMPLE')
    
    DX.Uice_SD.attrs = {'units':'m s-1', 
    'long_name':'Ensemble standard deviation of eastward sea ice drift velocity', } 
    DX.Vice_SD.attrs = {'units':'m s-1', 
    'long_name':'Ensemble standard deviation of northward sea ice drift velocity', }

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
    uvpairs = [('uice', 'vice'), ('uocean', 'vocean'), 
                ('Uice', 'Vice'), ('Uocean', 'Vocean')]

    uvstrs = ''

    for uvpair in uvpairs:

        if hasattr(DX, uvpair[0]) and hasattr(DX, uvpair[1]):
            uvc0_ = DX[uvpair[0]] + 1j* DX[uvpair[1]]
            uvc_rot = uvc0_ * np.exp(-1j*magdec_rad)

            DX[uvpair[0]].data = uvc_rot.real
            DX[uvpair[1]].data = uvc_rot.imag

            if hasattr( DX[uvpair[0]], 'heading_corr'):
                DX[uvpair[0]].attrs['heading_corr'] += '\n'+magdec_str
                DX[uvpair[1]].attrs['heading_corr'] += '\n'+magdec_str
            else:
                DX[uvpair[0]].attrs['heading_corr'] = magdec_str
                DX[uvpair[1]].attrs['heading_corr'] = magdec_str

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

    return DX

