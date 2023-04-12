""" Get the spectra """

import pandas as pd
import sys
import numpy as np
ddir = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data"
sys.path.append("/Users/annaho/Dropbox/astro/tools/Spectra")
from normalize import smooth_spec


def load_smoothed_spec(x, y, ivar):
    smoothed = smooth_spec(x, y, ivar, 1)
    return smoothed


def load_spec_1():
    """ read in the spectrum """
    inputf = "%s/opt/LRIS/lris20220923_improved_redux.spec" %ddir
    tab = pd.read_fwf(inputf, skiprows=np.arange(147))
    wl = tab['# wavelen'].values # AA
    # erg/cm2/s/Ang, absolute calibration approximate
    flam = tab['flux'].values 
    eflam[flam<0] = 1E8
    flam[flam<0] = 0
    eflam = tab['flux_unc'].values
    return wl, flam, eflam


def load_spec_2():
    """ read in the spectrum """
    #inputf = "%s/opt/LRIS/ZTF22abftjko_20221006_Keck1_v1.ascii" %ddir
    inputf = "%s/opt/LRIS/lris20221006_adjust.spec" %ddir
    tab = pd.read_fwf(inputf, skiprows=np.arange(150))
    wl = tab['wavelen'].values # AA
    # erg/cm2/s/Ang, absolute calibration approximate
    flam = tab['flux'].values 
    eflam[flam<0] = 1E8
    flam[flam<0] = 0
    eflam = tab['flux_unc'].values
    return wl, flam, eflam
