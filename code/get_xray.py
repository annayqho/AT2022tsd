""" Get the X-ray data of AT2022tsd """
import pandas as pd
import numpy as np
from astropy.time import Time
import vals

ddir = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/xray"

def load_swift_luminosity():
    """ load the Swift XRT data """
    df = pd.read_table(
            ddir+"/AT2022tsd_XRT_binned_luminosity.qdp",delimiter=" ")

    # Dates
    mjd1 = df['!col1'].values/3600/24+59856.387276
    t = Time(mjd1, format='mjd')
    Ls = df['col4']
    lLs = np.abs(df['col6'].values)
    uLs = df['col5']

    return t, Ls, lLs, uLs


def load_swift_counts():
    """ load the binned counts """
    df = pd.read_table(
            ddir+"/AT2022tsd_XRT_binned.qdp", 
            names=['t','t0','t1','ct','lct','uct'],delimiter='\t')
