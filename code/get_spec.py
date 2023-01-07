""" Get the spectra """

import pandas as pd
import numpy as np
ddir = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data"


def load_spec_1():
    """ read in the spectrum """
    inputf = "%s/opt/LRIS/lris20220923_improved_redux.spec" %ddir
    tab = pd.read_fwf(inputf, skiprows=np.arange(147))
    wl = tab['# wavelen'].values # AA
    flam = tab['flux'].values # erg/cm2/s/Ang, absolute calibration approximate
    eflam = tab['flux_unc'].values
    return wl, flam, eflam


def load_spec_2():
    """ read in the spectrum """
    inputf = "%s/opt/LRIS/ZTF22abftjko_20221006_Keck1_v1.ascii" %ddir
    tab = pd.read_fwf(inputf, skiprows=np.arange(150))
    wl = tab['wavelen'].values # AA
    flam = tab['flux'].values # erg/cm2/s/Ang, absolute calibration approximate
    eflam = tab['flux_unc'].values
    return wl, flam, eflam
