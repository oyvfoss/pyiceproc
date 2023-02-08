'''
Functions for making open water correction.
'''

from sigpyproc.sig_calc import runningstat, daily_average
from sigpyproc import sig_append, sig_draft

from scipy.interpolate import interp1d
import numpy as np
from matplotlib import pyplot as plt

def get_OWSD(DX, method = 'LE'):
    '''
    Get the surface depth during open water periods only.

    Returns DataArray containing the same variable as the input - 
    but with ice entries masked.
    '''
    OWSD = DX['SURFACE_DEPTH_%s'%method].where(
                            DX.ICE_IN_SAMPLE_ANY==False)
    return OWSD


def get_LP_OWSD(OWSD, thr_reject_from_net_median = 0.15, 
        min_frac_daily = 0.025, run_window_days = 3, ):
    '''
    Compute an estimate of the long-time averaged surface depth in open water
    (open water surface depth, OWSD).
    
    In an ideal case, OWSD should be equal to zero.

    Steps:
    1. Reject instances where OWSD deviates from the OWSD deployment median by 
       more than *thr_reject_from_net_median* (meters, default = 0.15).
    2. Compute ensemble median values of the OWSD resulting from (1).
    3. Compute daily medians of the ensemble means in (2).
       Reject days where less than *min_frac_daily* (default = 0.025) of the 
       ensembles contain open-water samples.
    4. Linearly interpolate between missing daily values to get a continuous
       daily time series.
    5. Smoothe this daily time series with a running mean of window length 
       *run_window_days* (default=3).
    '''
    # 1. Compute initial median and reject values away from the median 
    #    by *thr_reject_from_netmedian* [m]
    OWSD_full_median = OWSD.median()
    OWSD_filt = OWSD.where(
        np.abs(OWSD-OWSD_full_median)<thr_reject_from_net_median)
    
    # 2. Compute ensemble medians
    OWSD_med = OWSD_filt.median(dim = 'SAMPLE')
    
    #fig, ax = plt.subplots(2, 1, sharex = True)
    #ax[0].plot(OWSD.TIME, OWSD_med)

    # 3. Compute daily medians ()
    Ad, td = daily_average(OWSD_med, OWSD.TIME, min_frac = min_frac_daily, 
                         axis = -1, function = 'median')
    
    # 4. Interpolate to continuous function (daily)
    Ad_interp = interp1d(td.data[~np.isnan(Ad)], Ad[~np.isnan(Ad)], 
        bounds_error = False)(td.data)

    # 5. Smooth with running mean
    RS = runningstat(Ad_interp, run_window_days)
    
    # Export filtered, ensemble median, daily averaged, smoothed daily OWSD.
    # Also daily time array (td+0.5) of teh midpoint of the daily estimates.
    return RS['mean'], td+0.5


def get_Beta_from_OWSD(DX, 
        thr_reject_from_net_median = 0.15, 
        min_frac_daily = 0.025, run_window_days = 3,):
    '''
    Estimate sound speed correction BETA.
    '''

    # Obtain (all) estimates of daily, smoothed OWSDs
    OWSD_full_LE = get_OWSD(DX, method = 'LE')
    OWSD_full_AST = get_OWSD(DX, method = 'AST')

    # Obtain estimates of daily, smoothed OWSDs
    OWSD_LP_LE, _ = get_LP_OWSD(OWSD_full_LE, 
        thr_reject_from_net_median = thr_reject_from_net_median,
        min_frac_daily = min_frac_daily)
    OWSD_LP_AST, td = get_LP_OWSD(OWSD_full_AST, 
        thr_reject_from_net_median = thr_reject_from_net_median,
        min_frac_daily = min_frac_daily)


    # Obtain daily, smoothed instrument depths     
    depth_med = DX.depth.median(dim = 'SAMPLE')
    depth_med_daily, _ = daily_average(depth_med, DX.TIME, td = td-0.5, 
                                axis = -1, function = 'median')
    RS_depth = runningstat(depth_med_daily, run_window_days)
    depth_lp = RS_depth['mean']

    # Obtain Beta (sound speed correction factors)
    BETA_LE = depth_lp/(depth_lp - OWSD_LP_LE)
    BETA_AST = depth_lp/(depth_lp - OWSD_LP_AST)

    DX = sig_append.add_to_sigdata(DX, BETA_LE, td, 
            'BETA_open_water_corr_LE')
    DX = sig_append.add_to_sigdata(DX, BETA_AST, td, 
            'BETA_open_water_corr_AST')

    # Append the open water estimates as well
    DX = sig_append.add_to_sigdata(DX, OWSD_LP_LE, td, 
            'OW_surface_before_correction_LE')
    DX = sig_append.add_to_sigdata(DX, OWSD_LP_AST, td, 
            'OW_surface_before_correction_AST')
    return DX


def compare_OW_correction(DX, show_plots = True):
    '''
    Note: Run this *after* running *get_Beta_from_OWSD* but *before*
    running *sig_draft.calculate_draft()* again.
    '''

    DX0 = DX.copy()
    DX2 = DX.copy()
    DX2 = sig_draft.calculate_draft(DX2)


    print('LE: Mean (median) offset: %.1f cm (%.1f cm)'%(
        DX.OW_surface_before_correction_LE.mean()*1e2, 
        np.nanmedian(DX.OW_surface_before_correction_LE)*1e2))

    print('AST: Mean (median) offset: %.1f cm (%.1f cm)'%(
        DX.OW_surface_before_correction_AST.mean()*1e2, 
        np.nanmedian(DX.OW_surface_before_correction_AST)*1e2))

    print('LE: Mean (median) dBETA: %.1f (%.1f)'%(
        (DX.BETA_open_water_corr_LE.mean()-1)*1e3, 
        (np.nanmedian(DX.BETA_open_water_corr_LE)-1)*1e3))

    print('AST: Mean (median) dBETA: %.1f (%.1f)'%(
        (DX.BETA_open_water_corr_AST-1).mean()*1e3, 
        np.nanmedian((DX.BETA_open_water_corr_AST)-1)*1e3))


    print('LE - MEAN SEA ICE DRAFT:\n'+
          'Before OW correction: %.2f m'%DX0.SEA_ICE_DRAFT_MEDIAN_LE.mean()+
          '\nAfter OW correction: %.2f m'%DX2.SEA_ICE_DRAFT_MEDIAN_LE.mean())

    print('AST - MEAN SEA ICE DRAFT:\n'+
          'Before OW correction: %.2f m'%DX0.SEA_ICE_DRAFT_MEDIAN_AST.mean()+
          '\nAfter OW correction: %.2f m'%DX2.SEA_ICE_DRAFT_MEDIAN_AST.mean())


    # Figures
    if show_plots:
        fig, ax = plt.subplots(2, 1, sharex = True)

        ax[0].plot_date(DX2.TIME, DX.OW_surface_before_correction_LE, '-', label = 'LE')
        ax[0].plot_date(DX2.TIME, DX.OW_surface_before_correction_AST, '-', label = 'AST')

        ax[1].plot_date(DX2.TIME, DX.BETA_open_water_corr_LE, '-', label = 'LE')
        ax[1].plot_date(DX2.TIME, DX.BETA_open_water_corr_AST, '-', label = 'AST')

        for axn in ax: 
            axn.legend()
            axn.grid()
        labfs = 9
        ax[0].set_ylabel('Estimated open water\nsurface depth [m]', 
            fontsize = labfs)
        ax[1].set_ylabel('BETA (OWSD correction factor)', 
            fontsize = labfs)

        fig, ax = plt.subplots(1, 2, sharex = True, sharey = True)
        ax[0].scatter(DX0.time_average, DX0.SURFACE_DEPTH_LE, marker='.', 
                        color = 'k', alpha = 0.05, s = 0.3, label = 'Uncorrected')
        ax[0].scatter(DX.time_average, DX2.SURFACE_DEPTH_LE, marker='.', 
                        color = 'r', alpha = 0.05, s = 0.3, label = 'Corrected')

        ax[1].scatter(DX0.TIME, DX0.SEA_ICE_DRAFT_MEDIAN_LE, marker='.', 
                        color = 'k', alpha = 0.05, s = 0.3, label = 'Uncorrected')
        ax[1].scatter(DX.TIME, DX2.SEA_ICE_DRAFT_MEDIAN_LE, marker='.', 
                        color = 'r', alpha = 0.05, s = 0.3, label = 'Corrected')
        ax[0].set_title('LE Surface depth (ALL)')
        ax[1].set_title('LE sea ice draft (ice only, ensemble averaged)')

        for axn in ax: 
            axn.legend()
            axn.grid()
            axn.set_ylabel('[m]')

        labfs = 9
        ax[0].set_ylabel('Estimated open water\nsurface depth [m]', 
            fontsize = labfs)
        ax[1].set_ylabel('BETA (OWSD correction factor)', 
            fontsize = labfs)

        # Dummy for date axis..
        ax[0].plot_date(DX.time_average[0, 0],DX2.SURFACE_DEPTH_LE[0, 0] )
        
        ax[0].invert_yaxis()
        plt.show()

if False:
    OWSD_LP_LE, td = get_long_term_OWdist(OWDEP_LE, )
    print('LE: Mean (median) offset: %.1f cm (%.1f cm)'%(
        OWDEP_LP_LE.mean()*1e2, np.median(OWDEP_LP_LE)*1e2))

    OWDEP_LP_AST, td = get_long_term_OWdist(OWDEP_AST, )
    print('AST: Mean (median) offset: %.1f cm (%.1f cm)'%(
        OWDEP_LP_AST.mean()*1e2, np.median(OWDEP_LP_AST)*1e2))
