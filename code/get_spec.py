""" Get the spectra """

import pandas as pd
import sys
import numpy as np
ddir = "../data"
sys.path.append("/Users/annaho/Dropbox/astro/tools/Spectra")
from normalize import smooth_spec


def load_smoothed_spec(x, y, ivar):
    smoothed = smooth_spec(x, y, ivar, 1)
    return smoothed


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
