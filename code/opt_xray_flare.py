""" LC of the simultaneously observed optical and X-ray flare """

import pandas as pd
import numpy as np
from astropy.io import fits as pyfits
from get_opt import *
from get_xray import *

def flare_lc():
    """ Get the optical and X-ray light curves during the flare campaign """

    # The second LRIS flare was simultaneous with a Chandra observation
    start = 59942.42384
    end = 59942.43736

    # Get the optical points
    dat = get_full_opt()
    t = dat['mjdstart'].values
    choose = np.logical_and(t>=start, t<=end)
    dat_opt = dat[choose]

    # Get the relevant Chandra observation ID
    dat = load_chandra()
    obsid = 0
    for i,m in enumerate(dat['MJD'].values):
        if np.logical_and(m<=start, m+dat['Exp'].values[i]*1000/86400>=end):
            obsid = dat['OBSID'].values[i]

    # From the header of the repro_evt2.fits file, get the t0 in MJD and MET
    dd = "..//data/xray"
    head = pyfits.open(dd + "/%s/repro/acisf27643_repro_evt2.fits" %obsid)[0].header
    t0_mjd = Time(head['DATE-OBS'], format='isot').mjd

    # Get the light curve
    lc = np.loadtxt(dd + "/%s/repro/xray_flare_lc.txt" %obsid)
    dt = lc[:,0]
    ct = lc[:,1]
    ect = lc[:,2]

    # Convert time to MJD
    mjd = t0_mjd + dt/86400

    # Plot the optical fluxes in i-band
    choose = dat_opt['flt']=='i'
    xo = dat_opt['mjdstart'][choose].values+dat_opt['exp'][choose].values/2/86400
    exo = dat_opt['exp'][choose].values/2/86400
    yo = dat_opt['flux_extcorr'][choose].values
    eyo = dat_opt['unc_extcorr'][choose].values
    #plt.errorbar(xo, yo, xerr=exo, yerr=eyo, fmt='o', c='goldenrod')

    # Plot the X-ray count rate (500s bins)
    exx = 250/86400 # half the bin size
    xx = mjd+exx # center of the bin
    yx = ct
    eyx = ect
    #plt.errorbar(xx, yx, xerr=exx, yerr=eyx, fmt='o', c='k')

    return obsid,xo,exo,yo,eyo,xx,exx,yx,eyx
