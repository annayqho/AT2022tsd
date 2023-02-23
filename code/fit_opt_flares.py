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
        print("widths")
        print(popt[0], popt[3])
        print(np.sqrt(pcov[0,0]), np.sqrt(pcov[3,3]))

        print("offsets")
        print(popt[2], popt[5])
        print("difference")
        print((popt[2]-popt[5]))
        print("error on diff")
        print(np.sqrt(np.sqrt(pcov[2,2])**2+np.sqrt(pcov[5,5])**2))

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


if __name__=="__main__":
    # Fit the first LRIS flare
    out = get_full_opt()
    dat = out[np.logical_and(out['mjdstart']>59866.5, out['mjdstart']<59873)]
    x = (dat['mjdstart']-dat['mjdstart'].values[0])*24
    y = dat['flux']
    ey = dat['unc']

    # Get the g-band data
    tofit = np.logical_and(dat['flt']=='g',dat['#instrument']=='KeckI/LRIS') # exponential
    #tofit = choose # gaussian
    t0 = x[tofit].values[0]
    xgfit = x[tofit]-t0
    plt.errorbar(x[tofit]-t0, y[tofit], ey[tofit], fmt='o', c='green')
    ygfit = y[tofit]
    eygfit = ey[tofit]
    plt.errorbar(xgfit, ygfit, eygfit, fmt='o', c='green')

    tofit = np.logical_and(dat['flt']=='i',dat['#instrument']=='KeckI/LRIS') 
    xifit = x[tofit]-t0
    yifit = y[tofit]
    eyifit = ey[tofit]
    plt.errorbar(xifit, yifit, eyifit, fmt='o', c='goldenrod')

    # Try a dual-band exponential.
    xfit = np.hstack((xgfit, xifit)).astype(float)
    yfit = np.hstack((ygfit, yifit)).astype(float)
    eyfit = np.hstack((eygfit, eyifit)).astype(float)
    p0 = [2, 0.1, 2]
    popt, pcov = curve_fit(
            exp_gi, xfit, yfit, sigma=eyfit, absolute_sigma=True, p0=p0)
    xplt = np.linspace(min(xfit), max(xfit))
    yplt = exp_gi(np.hstack((xplt,xplt)), *popt)
    ygplt = yplt[0:int((len(yplt)+1)/2)]
    yiplt = yplt[int((len(yplt)+1)/2):]
    plt.plot(xplt, ygplt, c='green')
    plt.plot(xplt, yiplt, c='yellow')
    #print(popt, np.sqrt(pcov[1,1]))
    # Result: 3.6, 0.50, 1.66. clearly doesn't do a good job.

    # Try a Gaussian
    #p0 = [0.66, 4.44, 1.05] # for g-band
    #popt, pcov = curve_fit(
    #        gauss, xgfit, ygfit, sigma=eygfit, absolute_sigma=True, p0=p0, 
    #        maxfev=10000)
    #xplt = np.linspace(-1,2)
    #yg = gauss(xplt, *popt)
    #plt.plot(xplt, yg, c='green')
    #plt.plot(xplt, yg*2, c='yellow')
    #print(popt)
    #print(popt[1], np.sqrt(pcov[1,1]))


    plt.show()
