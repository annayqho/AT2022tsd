""" Fit functions to the various optical flares to try and characterize
their shapes """

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from get_opt import *
import vals

def gauss(x, sigma, A, b):
    """ Gaussian distribution """
    return A*np.exp(-(x-b)**2/(2*sigma**2))


def two_gaussians(x, sigma1, A1, b1, sigma2, A2, b2):
    """ Gaussian distribution """
    return A1*np.exp(-(x-b1)**2/(2*sigma1**2)) +\
           A2*np.exp(-(x-b2)**2/(2*sigma2**2))


def fred(x, a, b, c, d):
    """ Fast rise exponential decay """
    return a*np.exp(-b*(((x-c)/(d))+((d)/(x-c))))


def exp(x, A, tau):
    """ Fast rise exponential decay """
    return A*np.exp(-x/tau)


def exp_gi(x, A, tau, fac):
    """ exponential fit, for dual band (Keck) 
    same number of observations in both bands """
    xg = x[0:int((len(x)+1)/2)]
    yg = A*np.exp(-xg/tau)
    yi = yg*fac
    return np.hstack((yg, yi))


def pow_gi(x, A, tau, fac):
    """ power-law fit, for dual band (Keck) 
    same number of observations in both bands """
    xg = x[0:int((len(x)+1)/2)]
    yg = A*xg**(-tau)
    yi = yg*fac
    return np.hstack((yg, yi))


def gauss_gi(x, sigma, A, b, fac):
    xg = x[0:int((len(x)+1)/2)]
    yg = A*np.exp(-(xg-b)**2/(2*sigma**2))
    yi = yg*fac
    return np.hstack((yg, yi[1:]))


def gauss_const_gi(x, sigma, A, b, fac, offset):
    xg = x[0:int((len(x)+1)/2)]
    yg = offset+A*np.exp(-(xg-b)**2/(2*sigma**2))
    yi = yg*fac
    return np.hstack((yg, yi))


def gauss_exp_gi(x, sigma, A, b, fac, B, tau):
    xg = x[0:int((len(x)+1)/2)]
    yg = A*np.exp(-(xg-b)**2/(2*sigma**2)) + B*np.exp(-xg/tau)
    yi = yg*fac
    return np.hstack((yg, yi))


def gauss_exp(x, sigma, A, b, B, tau):
    """ Fit a Gaussian + exponential for a single band """
    return A*np.exp(-(x-b)**2/(2*sigma**2)) + B*np.exp(-x/tau)


def gauss_const(x, sigma, A, b, offset):
    """ Fit a Gaussian + constant for a single band """
    return A*np.exp(-(x-b)**2/(2*sigma**2)) + offset


def ultraspec_gband_small(ax, return_residuals=False):
    """ Fit the small ULTRASPEC g-band flare """
    out = get_full_opt()

    dat = out[np.logical_and(
        out['#instrument']=='TNT/ULTRASPEC', out['flt']=='g')]

    # Choose the first flare
    isflare = dat['isflare']
    t0 = np.min(dat['mjdstart'][isflare])
    xmin = dat['mjdstart']>t0-0.01
    xmax = dat['mjdstart']<t0+0.01
    tofit = np.logical_and(xmin, xmax)
    x = dat['mjdstart'][tofit]*24*60-dat['mjdstart'][tofit].values[0]*24*60
    y = dat['flux'][tofit]
    ey = dat['unc'][tofit]

    # Plot the flare
    ax.errorbar(x, y, ey, fmt='o', c=vals.gc, ms=4)

    # Fit
    p0 = [5, 20, 10]
    popt,pcov = curve_fit(gauss, x, y, sigma=ey, absolute_sigma=True, p0=p0)
    xplt = np.linspace(min(x), max(x), 1000)
    yplt = gauss(xplt, *popt)
    ax.plot(xplt,yplt,lw=1.5,c='k',zorder=10,
            label=r"$f_\nu \propto e^{-(x-b)^2/(2\sigma^2)}$")
    ax.text(0.95, 0.95, r"$\sigma=1.5\pm0.1$",ha='right',va='top',fontsize=8,
            transform=ax.transAxes)

    # Get the residuals
    ymod = gauss(x, *popt)

    if return_residuals:
        return x, y-ymod, ey


def ultraspec_gband_large(ax, show='exponential', return_residuals=False):
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
    keep = np.logical_and(x>t0-0.048, x<t0+0.1)
    nofit = np.logical_and(x>t0-0.04, x<t0-0.02)
    tofit = np.logical_and(keep, ~nofit)
    xfit = (x[tofit]-x[tofit].values[0])*24*60
    yfit = dat['flux'][tofit]
    eyfit = dat['unc'][tofit]
    xflare = (x[keep]-x[tofit].values[0])*24*60
    yflare = dat['flux'][keep]
    eyflare = dat['unc'][keep]
    
    # Fit the exponential
    p0 = [10, 1]
    popt,pcov = curve_fit(exp, xfit, yfit, sigma=eyfit, 
                          absolute_sigma=True, p0=p0)

    if show=='exponential':
        # Plot the exponential with its best fit
        plt.errorbar(xfit, yfit, eyfit, fmt='o', c=vals.gc, zorder=2, ms=4)
        plt.errorbar(
            xflare, yflare, eyflare, fmt='o', c='lightgrey', zorder=0, ms=4)
        xplt = np.linspace(min(xfit), max(xfit))
        yplt = exp(xplt, *popt)
        plt.plot(
                xplt,yplt,lw=1.5,c='k',zorder=5,
                label=r'$f_\nu\propto e^{-t/\tau}$')
        return popt
        # Best fit tau: 108 +/- 11 minutes

    elif show=='gaussian':
        # Divide out the exponential
        ymod = exp(xflare, *popt)
        p0 = [1, 2, 18, 1, 2, 30]
        popt, pcov = curve_fit(two_gaussians,
                xflare, yflare-ymod, sigma=eyflare/ymod, 
                absolute_sigma=True, p0=p0)
        ygausfit = yflare-ymod
        plt.errorbar(xflare, ygausfit, yerr=eyflare, fmt='o', 
                     c=vals.gc, zorder=2, ms=4)
        xplt = np.linspace(min(xflare), max(xflare), 1000)
        yplt = two_gaussians(xplt, *popt)
        ax.plot(
                xplt,yplt,lw=1.5,c='k', zorder=10, 
                label="Sum of two Gaussians")

        if return_residuals:
            ygausmod = two_gaussians(xflare, *popt)
            return xflare, ygausfit-ygausmod, eyflare


def ultraspec_gband():
    """ Do the analysis for the ULTRASPEC g-band light curve """
    # Initialize figure
    fig,axarr = plt.subplots(figsize=(7,7))

    # Plot the full g-band LC
    ax = plt.subplot(4,1,1)
    out = get_full_opt()
    dat = out[np.logical_and(
        out['#instrument']=='TNT/ULTRASPEC', out['flt']=='g')]
    x = dat['mjdstart']-dat['mjdstart'].values[0]
    y = dat['flux']
    ey = dat['unc']
    ax.errorbar(x*24,y,ey,fmt='o',c=vals.gc,ms=4)
    ax.axhline(y=0, lw=0.5, c='k', zorder=10)
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Flux (uJy)")

    # Plot the "noise" in the first 1.5 hours
    ax = plt.subplot(4,2,3)
    choose = x*24 < 1.5
    ax.hist(y[choose],color='grey')
    ax.set_ylabel("# Points")
    ax.set_xlabel(r"$f_\nu$")
    ax.text(0.02, 0.95, "Noise (first 1.5 hrs)", transform=ax.transAxes,va='top',fontsize=8)
    ax.axvline(x=0, color='k', lw=0.5)

    # Zoom in on the first blip and fit a Gaussian
    ax = plt.subplot(4,2,5)
    xres, yres, eyres = ultraspec_gband_small(
            ax, return_residuals=True) # Gaussian params
    ax.set_xlabel("Time (mins)")
    ax.set_ylabel("Flux (uJy)")
    ax.legend(fontsize=8, loc='upper left')

    # Plot the residuals of the fit
    ax = plt.subplot(4,2,7)
    ax.hist(yres, color='grey')
    ax.text(0.02, 0.95, "Fit Residuals", transform=ax.transAxes,va='top',fontsize=8)
    ax.axvline(x=0, color='k', lw=0.5)
    ax.set_xlabel(r"$f_\nu-f_{\nu,\mathrm{model}}$")
    ax.set_ylabel("# Points")

    # Zoom in on the second obvious flare and fit an exponential
    ax = plt.subplot(4,2,4)
    ultraspec_gband_large(ax) # exponential params
    ax.set_xlabel("Time (mins)")
    ax.set_ylabel("Flux (uJy)")
    ax.legend(fontsize=8, loc='upper right')
    ax.text(0.05,0.05,r'$\tau=110\pm11$ mins',transform=ax.transAxes, fontsize=8)
    #plt.axhline(y=1, lw=0.5, c='k')

    # Divide the exponential out and fit a double Gaussian
    ax = plt.subplot(4,2,6)
    xres, yres, eyres = ultraspec_gband_large(
            ax, show='gaussian', return_residuals=True)
    ax.legend(loc='upper right', fontsize=8)
    ax.text(
            0.02, 0.03, 
            '$\sigma_1=2.3$ mins, $\sigma_2=3.7$ mins, sep=$11.4$ mins',
            transform=ax.transAxes, fontsize=8)
    ax.set_ylim(-9,16)
    ax.set_xlabel("Time (mins)")
    ax.set_ylabel("Flux (uJy)")

    # Plot the residuals of the fit
    ax = plt.subplot(4,2,8)
    ax.hist(yres, color='grey')
    ax.text(0.02, 0.95, "Fit Residuals", transform=ax.transAxes,va='top',fontsize=8)
    ax.axvline(x=0, color='k', lw=0.5)
    ax.set_xlabel(r"$f_\nu-f_{\nu,\mathrm{model}}$")
    ax.set_ylabel("# Points")

    # Show
    plt.tight_layout()
    plt.savefig(
            "ultraspec_flare_fits.png", dpi=200, 
            bbox_inches='tight', pad_inches=0.1)
    plt.close()
    #plt.show()


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
    popt, pcov = curve_fit(
            gauss, np.array(xgfit).astype(float), 
            np.array(ygfit).astype(float), 
            sigma=np.array(eygfit).astype(float),
            absolute_sigma=True, p0=[40, 4.4, -31], maxfev=10000)
    xplt = np.linspace(min(xgfit[1:]), max(xgfit))
    yplt = gauss(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.gc)
    # i-band
    xfit_temp = np.append(xifit,-93.6)
    yfit_temp = np.append(yifit,-7.0)
    eyfit_temp = np.append(eyifit,4.0)
    popt, pcov = curve_fit(
            gauss, np.array(xfit_temp).astype(float), 
            np.array(yfit_temp).astype(float), 
            sigma=np.array(eyfit_temp).astype(float),
            absolute_sigma=True, p0=[40, 6.4, -31], maxfev=10000)
    xplt = np.linspace(min(xifit), max(xifit))
    yplt = gauss(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.ic)

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

    p0 = [0.1, 20]
    popt, pcov = curve_fit(
            exp, xifit, yifit, sigma=eyifit, 
            absolute_sigma=True, p0=p0, maxfev=10000)
    xplt = np.linspace(min(xifit), max(xifit))
    yplt = exp(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.ic)

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
    print(popt)

    # Gaussian+offset, varying color
    ax = axarr[2,0]
    ax.text(0.95, 0.95, 'Gaussian+offset, varying color', ha='right', va='top',
            transform=ax.transAxes, fontsize=8)
    p0 = [9, 2, 100, 1.7] 
    popt, pcov = curve_fit(
            gauss_const, xgfit[1:], ygfit[1:], sigma=eygfit[1:], 
            absolute_sigma=True, p0=p0, maxfev=10000)
    xplt = np.linspace(min(xgfit[1:]), max(xgfit))
    yplt = gauss_const(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.gc)

    p0 = [9, 2, 100, 1.7] 
    popt, pcov = curve_fit(
            gauss_const, xifit, yifit, sigma=eyifit, 
            absolute_sigma=True, p0=p0, maxfev=10000)
    xplt = np.linspace(min(xifit), max(xifit))
    yplt = gauss_const(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.ic)
    ax.set_xticklabels([0,5,10,15,20,25]) # shift

    for ax in axarr[:,0]:
        ax.set_ylabel(r"$f_\nu$ (uJy)")
    for ax in axarr[-1,:]:
        ax.set_xlabel("Minutes")

    plt.tight_layout()
    plt.savefig("lris_gi_fit.png", dpi=200)
    plt.close()


def imacs():
    """ Fit various functions to the Magellan flares """
    fig,axarr = plt.subplots(1,3,figsize=(7,2.5),sharex=True,sharey=True)

    # Get the flare from this day
    out = get_full_opt()
    dat = out[out['#instrument']=='Magellan/IMACS']
    x = (dat['mjdstart']-dat['mjdstart'].values[0])*24*60
    y = dat['flux']
    ey = dat['unc']

    for ax in axarr.flatten():
        ax.errorbar(x, y, ey, fmt='s', c=vals.gc, ms=4)
        #ax.set_yscale('log')
        #ax.set_yticks([2,3,4,6])
        #ax.set_yticklabels([2,3,4,6])

    # Top left: Gaussian fit to the bands separately
    ax = axarr[0]
    ax.text(0.98, 0.98, 'Gaussian', ha='right', va='top',
            transform=ax.transAxes, fontsize=8)
    # g-band
    p0 = [40, 4.4, -31]
    popt, pcov = curve_fit(gauss, x, y, sigma=ey, 
            absolute_sigma=True, p0=p0, maxfev=10000)
    xplt = np.linspace(min(x), max(x))
    yplt = gauss(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.gc)
    chisq = sum((y-gauss(x, *popt))**2/ey**2)
    dof = len(x)-len(p0)
    red_chisq = chisq/dof
    ax.text(0.98, 0.92, "$\chi^2$=%s" %np.round(red_chisq,1), 
            ha='right', va='top', transform=ax.transAxes, fontsize=8)

    # Next row: fast rise exponential decay
    ax = axarr[1]
    ax.text(0.98, 0.98, 'FRED', ha='right', va='top',
            transform=ax.transAxes, fontsize=8)
    p0 = [2E22, 20, -3, 40]
    popt, pcov = curve_fit(
            fred, x, y, sigma=ey, 
            absolute_sigma=True, p0=p0, maxfev=10000)
    xplt = np.linspace(min(x), max(x), 1000)
    yplt = fred(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.gc)
    chisq = sum((y-fred(x, *popt))**2/ey**2)
    dof = len(x)-len(p0)
    red_chisq = chisq/dof
    ax.text(0.98, 0.92, "$\chi^2$=%s" %np.round(red_chisq,1), 
            ha='right', va='top', transform=ax.transAxes, fontsize=8)

    # Gaussian+offset, varying color
    ax = axarr[2]
    ax.text(0.98, 0.98, 'Gaussian+offset', ha='right', va='top',
            transform=ax.transAxes, fontsize=8)
    p0 = [6, 50, 7, 1.7] 
    popt, pcov = curve_fit(
            gauss_const, x, y, sigma=ey, 
            absolute_sigma=True, p0=p0, maxfev=10000)
    xplt = np.linspace(min(x), max(x))
    yplt = gauss_const(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.gc)
    chisq = sum((y-gauss_const(x, *popt))**2/ey**2)
    dof = len(x)-len(p0)
    red_chisq = chisq/dof
    ax.text(0.98, 0.92, "$\chi^2$=%s" %np.round(red_chisq,1), 
            ha='right', va='top', transform=ax.transAxes, fontsize=8)

    axarr[0].set_ylabel(r"$f_\nu$ (uJy)")
    for ax in axarr:
        ax.set_xlabel("Minutes")

    plt.tight_layout()
    #plt.show()

    # The outcome is that the FRED clearly does the best job.
    plt.savefig("imacs_fit.png", dpi=200)
    plt.close()


if __name__=="__main__":
    """ Fit various functions to the LT flare """
    fig,axarr = plt.subplots(1,3,figsize=(7,2.5),sharex=True,sharey=True)

    # Get the flare from this day
    out = get_full_opt()
    dat = out[np.logical_and(out['mjdstart']>59929, out['mjdstart']<59929.9)]
    x = (dat['mjdstart']-dat['mjdstart'].values[0])*24*60
    y = dat['flux']
    ey = dat['unc']

    # Just pick the first few points
    x = x[0:4]
    y = y[0:4]
    ey = ey[0:4]

    for ax in axarr.flatten():
        ax.errorbar(x, y, ey, fmt='s', c=vals.gc, ms=4)

    # Top left: Gaussian fit 
    ax = axarr[0]
    ax.text(0.98, 0.98, 'Gaussian', ha='right', va='top',
            transform=ax.transAxes, fontsize=8)
    # g-band
    p0 = [40, 4.4, -31]
    popt, pcov = curve_fit(gauss, x, y, sigma=ey, 
            absolute_sigma=True, p0=p0, maxfev=10000)
    xplt = np.linspace(min(x), max(x))
    yplt = gauss(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.gc)
    chisq = sum((y-gauss(x, *popt))**2/ey**2)
    dof = len(x)-len(p0)
    red_chisq = chisq/dof
    ax.text(0.98, 0.92, "$\chi^2$=%s" %np.round(red_chisq,1), 
            ha='right', va='top', transform=ax.transAxes, fontsize=8)
    print(popt)

    # Next row: fast rise exponential decay
    ax = axarr[1]
    ax.text(0.98, 0.98, 'FRED', ha='right', va='top',
            transform=ax.transAxes, fontsize=8)
    p0 = [100, 1, 1, 1]
    popt, pcov = curve_fit(
            fred, x, y, sigma=ey, 
            absolute_sigma=True, p0=p0, maxfev=10000)
    xplt = np.linspace(min(x), max(x), 1000)
    yplt = fred(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.gc)
    #chisq = sum((y-fred(x, *popt))**2/ey**2)
    #dof = len(x)-len(p0)
    #red_chisq = chisq/dof
    #ax.text(0.98, 0.92, "$\chi^2$=%s" %np.round(red_chisq,1), 
    #        ha='right', va='top', transform=ax.transAxes, fontsize=8)

    # Gaussian+offset, varying color
    ax = axarr[2]
    ax.text(0.98, 0.98, 'Gaussian+offset', ha='right', va='top',
            transform=ax.transAxes, fontsize=8)
    p0 = [4, 11, 3, 1]
    popt, pcov = curve_fit(
            gauss_const, x, y, sigma=ey, 
            absolute_sigma=True, p0=p0, maxfev=10000)
    xplt = np.linspace(min(x), max(x))
    yplt = gauss_const(xplt, *popt)
    ax.plot(xplt, yplt, c=vals.gc)
    #chisq = sum((y-gauss_const(x, *popt))**2/ey**2)
    #dof = len(x)-len(p0)
    #red_chisq = chisq/dof
    #ax.text(0.98, 0.92, "$\chi^2$=%s" %np.round(red_chisq,1), 
    #        ha='right', va='top', transform=ax.transAxes, fontsize=8)

    axarr[0].set_ylabel(r"$f_\nu$ (uJy)")
    for ax in axarr:
        ax.set_xlabel("Minutes")

    plt.tight_layout()
    plt.show()

    # The outcome is that the FRED clearly does the best job.
    #plt.savefig("lt_fit.png", dpi=200)
    #plt.close()

