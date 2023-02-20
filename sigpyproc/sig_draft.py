'''
SIG_DRAFT.PY

Functions for calculating sea ice draft 
'''
import numpy as np

def calculate_draft(DX, corr_sound_speed_CTD = True, qual_thr = 8000,
                    LE_correction = 'AST'):
    '''
    Calculate ice draft.  


    If LE_correction = 'AST', the open water sound speed correction (if 
    available) of the LE-derived draft will be based on the AST open water 
    offset.   
    '''


    DX = calculate_surface_position(DX, qual_thr = qual_thr,
        corr_sound_speed_CTD = corr_sound_speed_CTD, le_ast ='LE')

    # Get surface position (LE and AST)
    for le_ast in ['LE', 'AST']:
        DX = calculate_surface_position(DX, qual_thr = qual_thr,
             corr_sound_speed_CTD = corr_sound_speed_CTD, le_ast =le_ast)

    # Reject LE measurements where LE diverges from AST by >0.5 m
    DX['SURFACE_DEPTH_LE_UNCAPPED'] = DX['SURFACE_DEPTH_LE']

    condition = np.abs(DX.SURFACE_DEPTH_LE - DX.SURFACE_DEPTH_AST)<0.5
    DX['SURFACE_DEPTH_LE'] = DX['SURFACE_DEPTH_LE'].where(condition)

    # Get sea ice draft (based on ice(no ice criteria))
    # Calculate ensemble medians
    for le_ast in ['LE', 'AST']:

        si_draft_ = DX['SURFACE_DEPTH_%s'%le_ast].data.copy()
        si_draft_[~DX.ICE_IN_SAMPLE.data] = np.nan

        DX['SEA_ICE_DRAFT_%s'%le_ast] = (('TIME', 'SAMPLE'), 
            si_draft_, {'long_name':'Sea ice draft at each sample (%s)'%le_ast, 
            'units':'m', 'note':DX['SURFACE_DEPTH_%s'%le_ast].note + 
                '\n\nSet to NaN where ICE_IN_SAMPLE==False'})

        DX['SEA_ICE_DRAFT_%s'%le_ast] = DX['SEA_ICE_DRAFT_%s'%le_ast].where(
                    DX['SEA_ICE_DRAFT_%s'%le_ast]>-0.3)

        DX['SEA_ICE_DRAFT_MEDIAN_%s'%le_ast] = (('TIME'), 
            np.nanmedian(si_draft_, axis=1),  {'long_name':
            'Median sea ice draft of each ensemble (%s)'%le_ast, 'units':'m', 
            'note':DX['SURFACE_DEPTH_%s'%le_ast].note + 
                '\n\nOnly counting instances with sea ice presence.'})

        DX['SEA_ICE_DRAFT_MEDIAN_%s'%le_ast] = DX['SEA_ICE_DRAFT_MEDIAN_%s'%le_ast].where(
                    DX['SEA_ICE_DRAFT_MEDIAN_%s'%le_ast]>-0.3)

    return DX



def calculate_surface_position(DX, corr_sound_speed_CTD = True, 
                                 qual_thr = 8000, le_ast ='AST', 
                                 LE_correction = 'AST'):
    '''
    Calculate distance between the surface measured by the altimeter 
    and the (mean) ocean surface.  

    If LE_correction = 'AST', the open water sound speed correction (if 
    available) of the LE-derived draft will be based on the AST open water 
    offset.   
    '''
    le_ast = le_ast.upper()
    if le_ast == 'AST':
        alt_dist_attr = 'Average_AltimeterDistanceAST'
        alt_qual_attr = 'Average_AltimeterQualityAST'
    else:
        alt_dist_attr = 'Average_AltimeterDistanceLE'
        alt_qual_attr = 'Average_AltimeterQualityLE'

    note_str = ('From %s altimeter distances.'
                '\n\nComputed with the function '
                'sig_draft.calculate_surface_position().'%le_ast)

    if hasattr(DX, 'sound_speed_CTD') and corr_sound_speed_CTD:
        # Ratio between observed and nominal sound speed
        sound_speed_ratio_obs_nom = (DX.sound_speed_CTD.data[:, np.newaxis]
            /DX.Average_Soundspeed.data)
        note_str += ('\n- Altimeter length recomputed using updated '
        'sound speed (*sound_speed_CTD* field)')
    else:
        sound_speed_ratio_obs_nom = 1
        note_str += ('\n- **NOT** USING OBSERVED SOUND SPEED\n'
           '  (*sound_speed* undefined) -> USING PRE-PROGRAMMED, NOMINAL VALUE!')    

    beta_key = 'BETA_open_water_corr_%s'%le_ast

    if LE_correction:
        beta_key = 'BETA_open_water_corr_AST'

    if hasattr(DX, beta_key):
        BETA_ow = (DX[beta_key].data[:, np.newaxis]
                     * np.ones(DX.depth.shape))
        note_str += ('\n- Using the open water correction factor to '
            'sound speed (*%s* field)'%beta_key)
    else:
        BETA_ow = 1
        note_str += ('\n- OPEN WATER CORRECTION **NOT** APPLIED!')    

    SURF_POS = DX.depth - (DX[alt_dist_attr]
                           * np.cos(np.pi*DX.tilt_Average/180)
                           * sound_speed_ratio_obs_nom
                           * BETA_ow)

    # APPLY QUALITY CRITERION
    SURF_POS = SURF_POS.where(DX[alt_qual_attr]>qual_thr)
    note_str += ('\n- Samples where %s>%i'%(alt_qual_attr, qual_thr) +
        ' were discarded.')

    # STORE AS VARIABLE
    DX['SURFACE_DEPTH_%s'%le_ast] = (('TIME', 'SAMPLE'), SURF_POS.data, {
        'long_name':'Depth of the scattering surface observed by'
        ' the Altimeter (%s)'%le_ast, 
        'units':'m', 'note':note_str})

    return DX


