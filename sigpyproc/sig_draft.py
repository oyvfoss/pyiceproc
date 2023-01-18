'''
SIG_DRAFT.PY

Functions for calculating sea ice draft 
'''
import numpy as np

def calculate_draft(DX):
    '''
    Calculate ice draft.  
    '''

    DX = calculate_surface_position(DX)

    si_draft = DX['SURFACE_DEPTH'].data.copy()
    si_draft[~DX.ICE_IN_SAMPLE.data] = np.nan

    DX['SEA_ICE_DRAFT'] = (('TIME', 'ENSEMBLE'), si_draft, 
        {'long_name':'Sea ice draft at each sample', 'units':'m', 
        'note':DX.SURFACE_DEPTH.note + 
            '\n\nSet to NaN where ICE_IN_SAMPLE==False'})

    DX['SEA_ICE_DRAFT_MEDIAN'] = (('TIME'), np.nanmedian(si_draft, axis=1), 
        {'long_name':'Median sea ice draft of each ensemble', 'units':'m', 
        'note':DX.SURFACE_DEPTH.note + 
            '\n\nOnly counting instances with sea ice presence.'})

    return DX



def calculate_surface_position(DX):
    '''
    Calculate distance between the surface measured by the altimeter 
    and the (mean) ocean surface.  
    '''

    note_str = ('From AST altimeter distances.'
                '\n\nComputed with the function '
                'sig_draft.calculate_surface_position().')

    if hasattr(DX, 'sound_speed_CTD'):
        # Ratio between observed and nominal sound speed
        sound_speed_ratio_obs_nom = (DX.sound_speed_CTD.data[:, np.newaxis]
            /DX.Average_Soundspeed.data)
        note_str += ('\n- Altimeter length recomputed using updated'
        'sound speed (*sound_speed_CTD* field)')
    else:
        sound_speed_ratio_obs_nom = 1
        note_str += ('\n- **NOT** USING OBSERVED SOUND SPEED\n'
           '  (*sound_speed* undefined) -> USING PRE-PROGRAMMED, NOMINAL VALUE!')    

    if hasattr(DX, 'c_s_corr_ow'):
        c_s_corr_ow = DX.c_s_corr_ow.data
        note_str += ('\n- Using the open water correction factor to '
            'sound speed (*cs_corr_ow* field)')
    else:
        c_s_corr_ow = 1
        note_str += ('\n- OPEN WATER CORRECTION **NOT** APPLIED!')    

    SURF_POS = DX.depth - (DX.Average_AltimeterDistanceAST
                           * np.cos(np.pi*DX.tilt_Average/180)
                           * sound_speed_ratio_obs_nom
                           * c_s_corr_ow)

    DX['SURFACE_DEPTH'] = (('TIME', 'ENSEMBLE'), SURF_POS.data, {
        'long_name':'Depth of the scattering surface observed by the Altimeter', 
        'units':'m', 'note':note_str})

    return DX