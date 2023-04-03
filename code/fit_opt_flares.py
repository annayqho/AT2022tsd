""" Fit functions to the various optical flares to try and characterize
their shapes """

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from get_opt import *
import vals


def ultraspec_gband_small():
    """ Get the small ULTRASPEC g-band flare """
    out = get_full_opt()

    dat = out[np.logical_and(
        out['#instrument']=='TNT/ULTRASPEC', out['flt']=='g')]

    # Choose the first flare
    isflare = dat['isflare']
    t0 = np.min(dat['mjdstart'][isflare])
    xmin = dat['mjdstart']>t0-0.006
    xmax = dat['mjdstart']<t0+0.006
    tofit = np.logical_and(xmin, xmax)
    x = dat['mjdstart'][tofit]*24*60-dat['mjdstart'][tofit].values[0]*24*60
    y = dat['flux'][tofit]
    ey = dat['unc'][tofit]

    return x.values, y.values, ey.values


def ultraspec_gband_large():
    """ Fit the large ULTRASPEC g-band flare """
    out = get_full_opt()

    # Get full LC
    dat = out[np.logical_and(
        out['#instrument']=='TNT/ULTRASPEC', out['flt']=='g')]
    x = dat['mjdstart']
    y = dat['flux']
    ey = dat['unc']

    # Get flare region
    isflare = dat['isflare']
    t0 = np.max(dat['mjdstart'][isflare])
    keep = np.logical_and(x>t0-0.05, x<t0+0.5)
    xflare = (x[keep]-x[keep].values[0])*24*60
    yflare = dat['flux'][keep]
    eyflare = dat['unc'][keep]

    return xflare.values, yflare.values, eyflare.values
    

def ultraspec_rband():
    """ get the ULTRASPEC r-band flare """
    out = get_full_opt()

    # Get full LC
    dat = out[np.logical_and(
        out['#instrument']=='TNT/ULTRASPEC', out['flt']=='r')]
    x = dat['mjdstart']
    y = dat['flux']
    ey = dat['unc']

    # Get flare region
    isflare = dat['isflare']
    t0 = np.max(dat['mjdstart'][isflare])
    keep = np.logical_and(x>t0-0.016, x<t0+0.001)
    xflare = (x[keep]-x[keep].values[0])*24*60
    yflare = dat['flux'][keep]
    eyflare = dat['unc'][keep]

    return xflare.values, yflare.values, eyflare.values


def lris_gi():
    """ Fit various functions to the LRIS g+i flares """
    fig,axarr = plt.subplots(3,2,figsize=(5,6),sharex=True,sharey=True)

    # Get the flare from this day
    out = get_full_opt()
    dat = out[np.logical_and(out['mjdstart']>59866.5, out['mjdstart']<59873)]
    x = (dat['mjdstart']-dat['mjdstart'].values[0])*24*60
    y = dat['flux']
    ey = dat['unc']

    # Get the g-band data
    #tofit = np.logical_and(dat['flt']=='g',dat['#instrument']=='KeckI/LRIS') # exponential
    tofit = dat['flt']=='g'
    t0 = x[tofit].values[0]
    xgfit = x[tofit]-t0
    ygfit = y[tofit]
    eygfit = ey[tofit]

    # Get the i-band data
    #tofit = np.logical_and(dat['flt']=='i',dat['#instrument']=='KeckI/LRIS')  # exponential
    tofit = dat['flt']=='i' 
    xifit = x[tofit]-t0
    yifit = y[tofit]
    eyifit = ey[tofit]

    for ax in axarr.flatten():
        ax.errorbar(xifit, yifit, eyifit, fmt='D', c=vals.ic, ms=4)
        ax.errorbar(xgfit[1:], ygfit[1:], eygfit[1:], fmt='s', c=vals.gc, ms=4)
        ax.set_yscale('log')
        ax.set_yticks([2,3,4,6])
        ax.set_yticklabels([2,3,4,6])

    # Top left: Gaussian fit to the bands separately
    ax = axarr[0,0]
    ax.text(0.95, 0.95, 'Gaussian, varying color', ha='right', va='top',
            transform=ax.transAxes, fontsize=8)
    # g-band
    p0 = [40, 4.4, -31]
    popt, pcov = curve_fit(
            gauss, np.array(xgfit).astype(float), 
            np.array(ygfit).astype(float), 
            sigma=np.array(eygfit).astype(float),
            absolute_sigma=True, p0=p0, maxfev=10000)
    xplt = np.linspace(min(xgfit[1:]), max(xgfit))
    yplt = gauss(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.gc)
    chisq = sum((ygfit[1:]-gauss(xgfit[1:], *popt))**2/eygfit[1:]**2)
    dof = len(xgfit[1:])-len(p0)
    red_chisq = chisq/dof
    ax.text(0.98, 0.8, "g-band, $\chi^2$=%s" %np.round(red_chisq,1), 
            ha='right', va='top', transform=ax.transAxes, fontsize=8)

    # i-band
    xfit_temp = np.append(xifit,-93.6)
    yfit_temp = np.append(yifit,-7.0)
    eyfit_temp = np.append(eyifit,4.0)
    p0 = [40, 6.4, -31]
    popt, pcov = curve_fit(
            gauss, np.array(xfit_temp).astype(float), 
            np.array(yfit_temp).astype(float), 
            sigma=np.array(eyfit_temp).astype(float),
            absolute_sigma=True, p0=p0, maxfev=10000)
    xplt = np.linspace(min(xifit), max(xifit))
    yplt = gauss(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.ic)

    chisq = sum((yifit-gauss(xifit, *popt))**2/eyifit**2)
    dof = len(xifit)-len(p0)
    red_chisq = chisq/dof
    ax.text(0.98, 0.92, "i-band, $\chi^2$=%s" %np.round(red_chisq,1), 
            ha='right', va='top', transform=ax.transAxes, fontsize=8)

    # Top right: Gaussian fit to the bands together
    ax = axarr[0,1]
    ax.text(0.95, 0.95, 'Gaussian, constant color', ha='right', va='top',
            transform=ax.transAxes, fontsize=8)
    xfit = np.hstack((xgfit, xifit)).astype(float)
    yfit = np.hstack((ygfit, yifit)).astype(float)
    eyfit = np.hstack((eygfit, eyifit)).astype(float)
    p0 = [40, 6.4, -31, 2]
    popt, pcov = curve_fit(
            gauss_gi, xfit, yfit, sigma=eyfit, absolute_sigma=True, 
            p0=p0, maxfev=10000)
    #popt = p0
    xplt = np.linspace(min(xfit[1:]), max(xfit))
    yplt = gauss_gi(np.hstack((xplt,xplt[1:])), *popt)
    ygplt = yplt[0:int((len(yplt)+1)/2)]
    yiplt = yplt[int((len(yplt)+1)/2):]
    ax.plot(xplt, ygplt, c=vals.gc)
    ax.plot(xplt[1:], yiplt, c=vals.ic)

    # Next row: exponential, single and dual 
    ax = axarr[1,0]
    ax.text(0.95, 0.95, 'Exponential, varying color', ha='right', va='top',
            transform=ax.transAxes, fontsize=8)
    p0 = [0.1, 20]
    popt, pcov = curve_fit(
            exp, xgfit[1:], ygfit[1:], sigma=eygfit[1:], 
            absolute_sigma=True, p0=p0, maxfev=10000)
    xplt = np.linspace(min(xgfit[1:]), max(xgfit))
    yplt = exp(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.gc)
    chisq = sum((ygfit[1:]-exp(xgfit[1:], *popt))**2/eygfit[1:]**2)
    dof = len(xgfit[1:])-len(p0)
    red_chisq = chisq/dof
    ax.text(0.98, 0.8, "g-band, $\chi^2$=%s" %np.round(red_chisq,1), 
            ha='right', va='top', transform=ax.transAxes, fontsize=8)

    p0 = [0.1, 20]
    popt, pcov = curve_fit(
            exp, xifit, yifit, sigma=eyifit, 
            absolute_sigma=True, p0=p0, maxfev=10000)
    xplt = np.linspace(min(xifit), max(xifit))
    yplt = exp(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.ic)
    chisq = sum((yifit-exp(xifit, *popt))**2/eyifit**2)
    dof = len(xifit)-len(p0)
    red_chisq = chisq/dof
    ax.text(0.98, 0.92, "i-band, $\chi^2$=%s" %np.round(red_chisq,1), 
            ha='right', va='top', transform=ax.transAxes, fontsize=8)

    # Try a dual-band exponential.
    ax = axarr[1,1]
    ax.text(0.95, 0.95, 'Exponential, constant color', ha='right', va='top',
            transform=ax.transAxes, fontsize=8)
    p0 = [0.1, 20, 2]
    popt, pcov = curve_fit(
            exp_gi, xfit[1:], yfit[1:], sigma=eyfit[1:], absolute_sigma=True, p0=p0)
    xplt = np.linspace(min(xfit[1:]), max(xfit))
    yplt = exp_gi(np.hstack((xplt,xplt)), *popt)
    ygplt = yplt[0:int((len(yplt)+1)/2)]
    yiplt = yplt[int((len(yplt)+1)/2):]
    ax.plot(xplt, ygplt, c=vals.gc)
    ax.plot(xplt, yiplt, c=vals.ic)

    # Try a dual-band Gaussian with an offset
    ax = axarr[2,1]
    ax.text(0.95, 0.95, 'Gaussian+offset, constant color', ha='right', va='top',
            transform=ax.transAxes, fontsize=8)
    p0 = [40, 4, 60, 2, 1] 
    popt, pcov = curve_fit(
            gauss_const_gi, xfit[1:], yfit[1:], sigma=eyfit[1:], 
            absolute_sigma=True, p0=p0, maxfev=10000)
    xplt = np.linspace(min(xfit[1:]), max(xfit))
    yplt = gauss_const_gi(np.hstack((xplt,xplt)), *popt)
    ygplt = yplt[0:int((len(yplt)+1)/2)]
    yiplt = yplt[int((len(yplt)+1)/2):]
    ax.plot(xplt, ygplt, c=vals.gc)
    ax.plot(xplt, yiplt, c=vals.ic)

    # Gaussian+offset, varying color
    ax = axarr[2,0]
    ax.text(0.95, 0.95, 'Gaussian+offset, varying color', ha='right', va='top',
            transform=ax.transAxes, fontsize=8)
    p0 = [9, 2, 100, 1.7] 
    print("g-band")
    popt, pcov = curve_fit(
            gauss_const, xgfit[1:], ygfit[1:], sigma=eygfit[1:], 
            absolute_sigma=True, p0=p0, maxfev=10000)
    print(popt[3], np.sqrt(pcov[3,3]))
    xplt = np.linspace(min(xgfit[1:]), max(xgfit))
    yplt = gauss_const(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.gc)
    chisq = sum((ygfit[1:]-gauss_const(xgfit[1:], *popt))**2/eygfit[1:]**2)
    dof = len(xgfit[1:])-len(p0)
    red_chisq = chisq/dof
    ax.text(0.98, 0.9, "g-band, $\chi^2$=%s" %np.round(red_chisq,1), 
            ha='right', va='top', transform=ax.transAxes, fontsize=8)

    print("i-band")
    p0 = [9, 2, 100, 1.7] 
    popt, pcov = curve_fit(
            gauss_const, xifit, yifit, sigma=eyifit, 
            absolute_sigma=True, p0=p0, maxfev=10000)
    print(popt[3], np.sqrt(pcov[3,3]))
    xplt = np.linspace(min(xifit), max(xifit))
    yplt = gauss_const(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.ic)
    ax.set_xticklabels([0,5,10,15,20,25]) # shift
    chisq = sum((yifit-gauss_const(xifit, *popt))**2/eyifit**2)
    dof = len(xifit)-len(p0)
    red_chisq = chisq/dof
    ax.text(0.98, 0.8, "i-band, $\chi^2$=%s" %np.round(red_chisq,1), 
            ha='right', va='top', transform=ax.transAxes, fontsize=8)

    for ax in axarr[:,0]:
        ax.set_ylabel(r"$f_\nu$ (uJy)")
    for ax in axarr[-1,:]:
        ax.set_xlabel("Minutes")

    plt.tight_layout()
    plt.show()
    #plt.savefig("lris_gi_fit.png", dpi=200)
    #plt.close()


def imacs():
    """ Fit various functions to the Magellan flares """
    # Get the flare from this day
    out = get_full_opt()
    dat = out[out['#instrument']=='Magellan/IMACS']
    x = (dat['mjdstart']-dat['mjdstart'].values[0])*24*60
    y = dat['flux_extcorr']
    ey = dat['unc_extcorr']
    return x.values, y.values, ey.values


def lt():
    """ Retrieve the extinction-corrected LT flare """
    out = get_full_opt()
    dat = out[np.logical_and(out['mjdstart']>59929, out['mjdstart']<59929.9)]
    x = (dat['mjdstart']-dat['mjdstart'].values[0])*24*60
    y = dat['flux_extcorr']
    ey = dat['unc_extcorr']
    filt = dat['flt']
    return filt.values, x.values, y.values, ey.values


def calc_t90(nu, x, y, ey):
    """ Calculate the T90 of a flare, aka the time in which 
    90% of the flux was recorded 

    Parameters
    ----------
    f: filter
    x: time values
    y: flux
    ey: uncertainty on flux
    """

    # First, get the peak luminosity
    ind = np.argmax(y)
    Lpeak = y[ind]*1E-6*1E-23*4*np.pi*vals.dL_cm**2 * nu

    # Interpolate onto a grid
    xvals = np.linspace(min(x), max(x), 1000)
    yvals = np.interp(xvals, x, y)

    # Get the cumulative fluence
    total = np.cumsum(yvals)[-1]

    # Fractional light curve as a function of time
    fvals = np.cumsum(yvals)/total

    # Identify where you cross 5% for the last time and 95% for the first time
    tstart= xvals[np.where(fvals<0.05)[0][-1]]
    tend = xvals[np.where(fvals>0.95)[0][0]]

    # Calculate T90
    t90 = tend-tstart
    print(t90, Lpeak)
    plt.plot(xvals, fvals, c='grey')
    plt.axvline(x=tstart, c='k')
    plt.axvline(x=tend, c='k')
    plt.show()

    return t90, Lpeak


if __name__=="__main__":
    #f, x, y, ey = lt()
    #x, y, ey = imacs()
    #x, y, ey = ultraspec_gband_large()
    x, y, ey = ultraspec_rband()
    nu = 3E18/vals.ztf_pivot['r']
    t90, lpeak = calc_t90(nu, x, y, ey) 
    # LT: 10.4 minutes, 3.6E43 erg/s
    # IMACS: 16.1 minutes, 1.5E44 erg/s
    # ULTRASPEC g-band, short: 7.7 minutes. 1.8E43 erg/s
    # ULTRASPEC g-band, large: 78.0 minutes, 2.7E43 erg/s
    # ULTRASPEC r-band: 19.1 minutes, 3.0E43 erg/s
