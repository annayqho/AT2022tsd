""" Get the spectra """

import pandas as pd
import sys
import numpy as np
from scipy.signal import savgol_filter
ddir = "../../data"


def bin_spec(x, y, ey, bin_size):
    x_binned = []
    y_binned = []
    ey_binned = []

    for ii,x_val in enumerate(x):
        choose = np.abs(x-x_val)<bin_size
        if sum(choose)==1:
            x_binned.append(x_val)
            y_binned.append(y[choose][0])
            ey_binned.append(ey[choose][0])
        elif sum(choose)>1:
            mean,wsum = np.average(
                y[choose], weights=1/ey[choose]**2, returned=True)
            efmean = np.sqrt(1/wsum)
            x_binned.append(np.average(x[choose]))
            y_binned.append(mean)
            ey_binned.append(efmean)

    x_binned = np.array(x_binned)
    y_binned = np.array(y_binned)
    ey_binned = np.array(ey_binned)

    return x_binned,y_binned,ey_binned


def load_binned_spec(x, y, ey):
    binned = bin_spec(x, y, ey, 3)
    return binned


def load_spec(inputf, ran):
    tab = pd.read_fwf(inputf, skiprows=np.arange(ran))
    wl = tab['# wavelen'].values # AA
    # erg/cm2/s/Ang, absolute calibration approximate
    flam = tab['flux'].values 
    eflam = tab['flux_unc'].values

    # Set bad values to 0, with large errors
    eflam[flam<0] = 1E8
    eflam[np.isnan(flam)] = 1E8
    flam[flam<0] = 0
    flam[np.isnan(flam)] = 0
    return wl, flam, eflam


def load_spec_1():
    """ read in the spectrum """
    inputf = "%s/opt/LRIS/spectroscopy/lris20220923_improved_redux.spec" %ddir
    wl, flam, eflam = load_spec(inputf, 147)
    return wl, flam, eflam


def load_spec_2():
    """ read in the spectrum """
    inputf = "%s/opt/LRIS/spectroscopy/lris20221006_adjust.spec" %ddir
    wl, flam, eflam = load_spec(inputf, 150)
    return wl, flam, eflam
