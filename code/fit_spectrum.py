""" Fit the Keck/LRIS spectrum of AT2022tsd 

Fit the Galactic extinction corrected LRIS spectrum
MILES library (FWHM = 2.5 AA)
Commonly observed galaxy emission lines
    Halpha
    Hbeta
    Hgamma
    OII
    SII
    OIII
    OI
    NII
    Doublets fixed at theoretical flux ratio of 3:
        OI (6300,6364), OIII (4959, 5007), NII (6548, 6583)
Calculate line ratios, uncertainties by performing 1E4 Monte Carlo trials
using measured flux uncertainties

Use the 10-06 spectrum since it has clearer lines. This is very close to the
NOT photometry, which was on 10-04. So, let's use that as a reference.
Those filters are sdssr, sdssi, and sdssg. 

The spectra we have still have a lot of transient flux, so there's no point
in trying to fit the continuum. Just solve for the redshift for now.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from astropy import units as u
from specutils import Spectrum1D
import extinction
from dustmaps.sfd import SFDQuery
import warnings
from specutils.fitting import fit_generic_continuum
from astropy.coordinates import SkyCoord
import vals 

# Direc with all the data we'll use in this code
ddir = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data"
c = 3e18 # speed of light in AA


def gaussian(x, amp, cen, wid):
    """ Fit a Gaussian """
    return amp * np.exp(-(x-cen)**2 / wid)


def fit_redshift(wl_all, flux_all, eflux_all, l0, z0, window):
    """
    Given a galaxy spectrum, get the best-fit redshift

    Parameters
    ----------
    wl_all: the wl array of the spectrum
    flux_all = the flux array of the spectrum
    l0: the rest wavelength of the line in question
    z0: the initial guess for the redshift

    Returns
    -------
    z: best-fit redshift
    zerr: uncertainty on best-fit redshift
    """
    # Get rid of any nans
    wl_all = wl_all[~np.isnan(flux_all)]
    eflux_all = eflux_all[~np.isnan(flux_all)]
    flux_all = flux_all[~np.isnan(flux_all)]

    # Extract the region of interest
    lm = l0 * (z0 + 1)

    choose = np.logical_and(
            wl_all > (lm - window),
            wl_all < (lm + window))
    wl = wl_all[choose]
    flux = flux_all[choose]
    eflux = eflux_all[choose]

    #plt.step(wl_all,flux_all,lw=0.5,c='grey')
    #plt.plot(wl, flux)
    #plt.axhline(y=np.median(flux))

    # local continuum
    cont_range = np.logical_or(
            np.logical_and(wl_all > (lm - window),
                            wl_all < (lm -2)),
            np.logical_and(wl_all > (lm + 2),
                            wl_all < (lm + window)))
    cont = np.median(flux_all[cont_range])
    #plt.axhline(y=cont)
    #plt.show()
    
    # Normalize
    flux_norm = (flux / cont)-1

    # Fit the Gaussian to this region
    init_vals = [4, lm, 1]
    best_vals, covar = curve_fit(
            gaussian, wl, flux_norm, p0=init_vals, 
            sigma=eflux, absolute_sigma=True)

    # Plot the fit
    # plt.axvline(x=lm, lw=0.5, c='grey') # guess
    # plt.plot(wl, flux_norm, c='k')
    # xfit = np.linspace(min(wl), max(wl), 1000)
    # yfit = gaussian(xfit, best_vals[0], best_vals[1], best_vals[2])
    # plt.plot(xfit, yfit, c='r')

    # Convert the best-fit into an actual redshift
    center = best_vals[1]
    ecenter = np.sqrt(covar[1][1])
    z = (center-l0)/l0
    ez = ecenter/l0

    #plt.axvline(x=l0 * (z + 1), c='red', lw=0.5)
    #plt.show()

    # Return values
    return z, ez


def get_rest_wl():
    """ Get the rest wavelengths in air of the relevant lines 
    http://astronomy.nmsu.edu/drewski/tableofemissionlines.html """
    wl_lines = {}
    wl_lines['ha'] = [6562.819]
    wl_lines['heii'] = [4686]
    wl_lines['hei'] = [5875]
    wl_lines['hb'] = [4861.333]
    wl_lines['oii'] = [3726.032, 3728.815]
    wl_lines['oiii'] = [4958.911, 5006.843]
    wl_lines['nii'] = [6548.050, 6583.460]
    wl_lines['sii'] = [6716.440, 6730.810]
    #hgamma = 4340.471
    return wl_lines


def gaussian_weight_matrix(wl, L):
    """ Matrix of Gaussian weights 

    Parameters
    ----------
    wl: numpy ndarray
        pixel wavelength values
    L: float
        width of Gaussian in pixels

    Return
    ------
    Weight matrix
    """
    return np.exp(-0.5*(wl[:,None]-wl[None,:])**2/L**2)


def smooth_spec(wl, flux, ivar, L):
    """ Smooth a spectrum with a running Gaussian

    Parameters
    ----------
    wl: numpy ndarray
        pixel wavelength values
    flux: numpy ndarray
        pixel flux values
    ivar: numpy ndarray
        pixel inverse variances
    L: width of Gaussian in Angstroms

    Returns
    ------
    smoothed: numpy ndarray
        smoothed flux values
    """
    w = gaussian_weight_matrix(wl, L)
    bot = np.dot(ivar, w.T)
    top = np.dot((flux*ivar), w.T)
    bad = bot == 0
    smoothed = np.zeros(top.shape)
    smoothed[~bad] = top[~bad] / bot[~bad]
    return smoothed


def load_spec():
    """ read in the spectrum """
    inputf = "%s/opt/LRIS/lris20221006_adjust.spec" %ddir
    tab = pd.read_fwf(inputf, skiprows=np.arange(150))
    wl = tab['wavelen'].values # AA
    flam = tab['flux'].values # erg/cm2/s/Ang, absolute calibration approximate
    eflam = tab['flux_unc'].values
    return wl, flam, eflam


def ext_corr(wl, flam):
    """ correct for Galactic extinction
    from the IRSA service, Schlafly & Finkbeiner """
    R_V = 3.1
    A_V = 0.7413
    ext = extinction.fitzpatrick99(wl, A_V, R_V)
    fnu = 3.34E4 * wl**2 * flam # Jy
    mAB = -2.5*np.log10(fnu)+8.90 
    mAB_corr = mAB-ext
    fnu_corr = 10**((mAB_corr-8.90)/(-2.5))
    flam_corr = fnu_corr / 3.34E4 / wl**2
    return wl, flam_corr


def get_mags():
    """ Get the measured magnitudes, no extinction """
    mag_ref = {}
    for filt in ['g', 'r', 'i']:
        phot = pd.read_csv("%s/opt/opt.txt" %ddir)
        is_ref = phot['Date']==59856.0
        phot = phot[is_ref]
        mag_ref[filt] = phot['Mag'][phot['Filt']==filt].values[0]
    return mag_ref


def get_spec_mags():
    """ Get the observed raw spec magnitudes """
    spec_obs_mAB = {}
    for filt in ['g', 'r', 'i']:
        trans = np.loadtxt("%s/SLOAN_SDSS.%s.dat" %(ddir,filt))
        R = np.interp(wl, trans[:,0], trans[:,1]) # interpolate transmission func
        w_flam = R*flam # erg/cm2/s/AA
        choose = w_flam > 0 # for integration
        plt.plot(wl[choose], w_flam[choose])
        top = np.trapz(wl[choose]*w_flam[choose], wl[choose])
        bottom = np.trapz(wl[choose]*R[choose], wl[choose])
        mean_flam = top / bottom # erg/cm2/s/AA
        print(mean_flam)
        mean_fnu = 3.34E4 * (vals.sdss_pivot[filt])**2 * mean_flam # Jy
        print(mean_fnu)
        spec_obs_mAB[filt] = -2.5*np.log10(mean_fnu)+8.90 
    return spec_obs_mAB


def run_smooth_spec():
    """ Apply the spectrum smoothing. Currently not working... """
    wl, flam, eflam = load_spec()
    eflam[flam < 0] = 1E8
    flam[flam < 0] = min(flam[flam>0])
    flam = flam / 1E-17
    eflam = eflam / 1E-17
    ivar = 1/eflam**2
    ivar[ivar<1] = 1

    smoothed_spec = smooth_spec(wl, flam, ivar, L=10)
    norm_flux = flam / smoothed_spec
    norm_ivar = ivar * smoothed_spec**2
    plt.plot(wl ,flam, c='grey', lw=0.5)
    plt.plot(wl, smoothed_spec, c='red', lw=0.5)


def fit_z():
    """ Fit for the redshift """

    # Load spectrum 
    wl,f,ef = load_spec()

    # Initial guess
    z0 = 0.2565

    # Window size, in angstroms
    window = 20

    # Fit for centroid of a line
    wl_lines = get_rest_wl()
    zall = []
    ezall = []
    for key,val in wl_lines.items():
        if key=='sii':
            use = val[0:-1] # only fit the first SII line
        elif key=='oii':
            use = [] # doublet, makes things confusing
        else:
            use = val
        for line in use:
            #plt.axvline(x=line*(1+z0), c='red', lw=0.5)
            z, ez = fit_redshift(wl, f, ef, line, z0, window)
            zall.append(z)
            ezall.append(ez)
            print(key, z, ez)
    zall = np.array(zall)
    ezall = np.array(ezall)
    w = 1/ezall**2
    zmean = np.average(zall, weights=w)
    ezmean = np.std(zall)
    print("%s +/- %s" %(np.round(zmean,7), np.round(ezmean, 7)))


if __name__=="__main__":
    fit_z()

