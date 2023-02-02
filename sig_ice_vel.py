'''
SIG_ICE_VEL.PY

Functions for processing sea ice draft.
'''

def rotate_vels_magdec(DX):
    '''
    Rotate ocean and ice velocities to account for magnetic declination.

    Only rotates processed fields (uocean, vocean, uice, vice) -
    not raw variables (Average_VelNorth, AverageIce_VelNorth, etc..).
    '''
    
    assert hasattr(DX, 'magdec'), ("Didn't find magnetic declination (no" 
        " *magdec* attribute). Run sig_append.append_magdec()..")

    # Convert to radians
    magdec_rad = DX.magdec * np.pi/180
    # Make a documentation string (to add as attribute)
    magdec_mean = DX.magdec.mean().data
    magdec_str = ('Rotated CW by an average of %.2f degrees'%magdec_mean
                    ' to correct for magntic declination.')

    # Loop through different (u, v) variable pairs and rotate them
      
    uvpairs = [('uice', 'vice'), ('uocean', 'vocean'), 
                ('Uice', 'Vice'), ('Uocean', 'Vocean')]

    for uvpair in uvpairs:
        if hasattr(DX, uvpair[0]) and hasattr(DX, uvpair[1]):
            uvc0_ = DX[uvpair[0]] + 1j* DX[uvpair[1]]
            uvc_rot = uvc0_ * np.exp(-1j*magdec_rad)

            DX[uvpair[0]] = uvc_rot.real
            DX[uvpair[0]].attrs['heading_corr'] = magdec_str
            DX[uvpair[1]] = uvc_rot.imag
            DX[uvpair[1]].attrs['heading_corr'] = magdec_str

    return DX

