""" Plot the Chandra X-ray flares """

import numpy as np
import pandas as pd
from astropy.io import fits as pyfits
from astropy.time import Time
import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
from get_xray import *


def plot_flares(axarr):
    dd = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/xray"

    df = load_chandra()
    oids = df['OBSID'].values.astype(str)
    t = Time(df['MJD'].values, format='mjd')
    exp = (df['Exp'].values*1000)/86400 # days
    dt_start = (t.jd-vals.t0)/(1+vals.z)
    dt_dur = exp/(1+vals.z)
    dates = dt_start + dt_dur/2 # halfway through

    for i,oid in enumerate(oids):
        ax = axarr.flatten()[i]
        ff = dd + "/" + oid + "/repro/xray_flare_lc.txt"
        bf = dd + "/" + oid + "/repro/xray_bkg_lc.txt"
        ff2 = dd + "/" + oid + "/repro/xray_src_offset_lc.txt"

        dat = np.loadtxt(ff)
        ax.errorbar(
                dat[:,0]/60, dat[:,1], yerr=dat[:,2], marker="s", c='k',
                label='AT2022tsd',zorder=2,ms=3, lw=0.5)

        # Scale the count rate by the ratio of the size of the background
        # to the size of the source
        src_size = 2
        with open(dd + "/" + oid + "/repro/bkg_offset.reg") as inputf:
            l = inputf.readlines()[-1]
            bkg_size = float(l.split(',')[-1][0:-2])/2

        dat = np.loadtxt(bf)
        ax.plot(dat[:,0]/60,dat[:,1]*(src_size**2/bkg_size**2),c='lightgrey',
                lw=2, label='Background', zorder=0)

        if i==0:
            ax.legend(loc='upper left', fontsize=8)

        dat = r"$\Delta t=%s$" %np.round(dates[i], 1)
        if np.logical_or(i==0, i==2):
            ax.text(
                    0.88, 0.95, dat, ha='right', va='top', 
                    fontsize=8, transform=ax.transAxes)
        else:
            ax.text(
                    0.98, 0.95, dat, ha='right', va='top', 
                    fontsize=8, transform=ax.transAxes)

    axarr.flatten()[-1].axis('off')

    for ax in axarr[:,0]:
        ax.set_ylabel("Count Rate")
    for ax in axarr[-1,:]:
        ax.set_xlabel("Minutes")
