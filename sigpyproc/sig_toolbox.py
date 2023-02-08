
import numpy as np
import matplotlib.pyplot as plt

def plot_ellipse_icevel(DX, lp_days = 5, ax = None, 
                return_ax = True):
    '''
    Plot of u and v ice drift components
    low pass filtered with a running mean of *lp_days*.

    Showing the mean current vector, the low-pass-filtered 
    and subsampled currents, and the semi-major and -minor axes
    of the variance ellipse. 
    '''
    assert hasattr(DX, 'Uice'), (
        'No "Uice" field. Run sig_ice_vel.calculate_drift()..')

    print('ELLIPSE PLOT: Interpolate over nans.. \r', end = '')

    uip = DX.Uice.interpolate_na(dim = 'TIME', limit=10).data
    vip = DX.Vice.interpolate_na(dim = 'TIME', limit=10).data

    print('ELLIPSE PLOT: Low pass filtering..    \r', end = '')
    # LPFed 
    wlen = int(np.round(lp_days/(DX.Plan_ProfileInterval/60/60/24)))
    ULP = np.convolve(uip, np.ones(wlen)/wlen,
                mode = 'valid')[::wlen]
    VLP = np.convolve(vip, np.ones(wlen)/wlen,
                mode = 'valid')[::wlen]

    print('ELLIPSE PLOT: Calculating ellipse (from LPed data).. \r', end = '')

    # Ellipse
    thp, majax, minax = _uv_angle(
        ULP - np.nanmean(ULP), VLP - np.nanmean(VLP))

    # Mean
    UM, VM = np.nanmean(ULP), np.nanmean(VLP)


    print('ELLIPSE PLOT: Plotting..                              \r', end = '')

    if ax == None:
        fig, ax = plt.subplots(figsize = (10, 10))

    ax.set_aspect('equal')

    ax.plot(uip, vip, '.', ms = 1, color = 'Grey', 
        alpha = 0.3, lw = 2, zorder = 0)
    ax.plot(DX.Uice.data[-1], DX.Vice.data[-1], '.', ms= 1,color = 'k', alpha =0.5,
            lw = 2, label ='Full')

    ax.plot(ULP, VLP, '.', ms = 3, color = 'b', alpha = 0.5)
    ax.plot(ULP[0], VLP[0], '.', ms = 5, color = 'b', alpha = 0.5, 
        label = '%.1f-day means'%(lp_days), zorder = 0)

    vmaj = np.array([-majax*np.sin(thp), majax*np.sin(thp)])
    umaj = np.array([-majax*np.cos(thp), majax*np.cos(thp)])
    vmin = np.array([-minax*np.sin(thp+np.pi/2),
                        minax*np.sin(thp+np.pi/2)])
    umin = np.array([-minax*np.cos(thp+np.pi/2),
                        minax*np.cos(thp+np.pi/2)])

    ax.plot(UM + umaj , VM + vmaj, '-k', lw = 2, label ='Maj axis')
    ax.plot(UM + umin , VM + vmin, '--k', lw = 2, label ='Min axis')

    ax.quiver(0, 0, UM, VM, 
        color = 'r', scale_units = 'xy', scale = 1, width = 0.03, 
        headlength = 2, headaxislength = 2, alpha = 0.6, 
        label = 'Mean (u: %.2f, v: %.2f)'%(UM, VM),  edgecolor='k', 
        linewidth = 0.6)

    ax.set_ylabel('v [m s$^{-1}$]'); ax.set_xlabel('u [m s$^{-1}$]')
    ax.legend(fontsize= 10, loc =3, handlelength = 1, ncol = 2) 

    ax.set_title('Ice drift velocity components')
    ax.grid()
    plt.show()

    if return_ax:
        return ax





def _uv_angle(u, v):
    '''
    Finds the principal angle angle in [-pi/2, pi/2] where the squares
    of the normal distance  to u,v are maximised. Ref Emery/Thompson
    pp 327.

    Also returns the standard deviation along the semimajor and
    semiminor axes.

    '''
    if np.nanmean(u)>1e-7 or np.nanmean(v)>1e-7:
        print('Mean of u and/or v is nonzero. Removing mean.')
        u = u - np.nanmean(u) ; v = v - np.nanmean(v)
    thp = 0.5* np.arctan2(2*np.nanmean(u*v),(np.nanmean(u**2)-\
                          np.nanmean(v**2))) #ET eq 4.3.23b
    uvcr = (u + 1j * v) * np.exp(-1j * thp)
    majax = np.nanstd(uvcr.real)
    minax = np.nanstd(uvcr.imag)

    return thp, majax, minax