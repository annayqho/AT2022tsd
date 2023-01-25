""" Get the X-ray data of AT2022tsd """
import pandas as pd
import numpy as np
from astropy.time import Time
import vals

ddir = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/xray"


def load_swift_counts():
    """ load the unbinned counts """
    df = pd.read_table(
            ddir+"/AT2022tsd_XRT_unbinned.qdp", delimiter='\t',
            skiprows=np.hstack((np.arange(0,14),np.array([25,26]))))
    return df


def load_swift():
    """ Load the full counts, flux, luminosity """
    df = load_swift_counts()

    # Convert counts to flux
    conv = 5.10E-11 # From the XRT spectrum fit
    df['Flux'] = df['Rate    ']*conv
    df['Fluxpos'] = df['Ratepos ']*conv
    df['Fluxneg'] = df['Rateneg']*conv

    # Now the luminosity
    df['L'] =  df['Flux']*4*np.pi*(vals.dL_cm)**2
    df['Lpos'] =  df['Fluxpos']*4*np.pi*(vals.dL_cm)**2
    df['Lneg'] =  df['Fluxneg']*4*np.pi*(vals.dL_cm)**2

    return df
