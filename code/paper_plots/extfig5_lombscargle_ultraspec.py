""" Do time-series analysis of the optical light-curve data """

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.timeseries import LombScargle
import sys
sys.path.append("..")
import vals
from get_opt import *
#from opt_lc import plot_det, plot_lim
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 7 # The maximum allowed for ED figures


def get_gband_flare():
    dat = get_full_opt()
    choose = np.logical_and(
            dat['#instrument']=='TNT/ULTRASPEC', dat['flt']=='g')
    dat = dat[choose]
    f = dat['flux'].values.astype(float)
    ef = dat['unc'].values.astype(float)
    dt = np.array(dat['mjdstart'].values-dat['mjdstart'].values[0]).astype(float)

    # Parameters
    padding = 0.01
    sig = 4

    isdet = f/ef>=sig
    flare_start = min(dt[isdet])-padding
    flare_end = max(dt[np.logical_and(isdet, dt<0.7)])+padding
    choose = np.logical_and(dt>=flare_start, dt<=flare_end)
    x = dt[choose]
    y = f[choose]
    ey = ef[choose]
    return x,y,ey


def get_gband_noise():
    dat = get_full_opt()
    choose = np.logical_and(
            dat['#instrument']=='TNT/ULTRASPEC', dat['flt']=='g')
    dat = dat[choose]
    f = dat['flux'].values.astype(float)
    ef = dat['unc'].values.astype(float)
    dt = np.array(dat['mjdstart'].values-dat['mjdstart'].values[0]).astype(float)
    choose = dt < 0.065
    x = dt[choose]
    y = f[choose]
    ey = ef[choose]
    return x,y,ey


def get_gband_all():
    dat = get_full_opt()
    choose = np.logical_and(
            dat['#instrument']=='TNT/ULTRASPEC', dat['flt']=='g')
    dat = dat[choose]
    f = dat['flux'].values.astype(float)
    ef = dat['unc'].values.astype(float)
    dt = np.array(dat['mjdstart'].values-dat['mjdstart'].values[0]).astype(float)
    return dt, f, ef


def get_rband_flare():
    dat = get_full_opt()
    choose = np.logical_and(
            dat['#instrument']=='TNT/ULTRASPEC', dat['flt']=='r')
    dat = dat[choose]
    f = dat['flux'].values.astype(float)
    ef = dat['unc'].values.astype(float)
    dt = np.array(dat['mjdstart'].values-dat['mjdstart'].values[0]).astype(float)

    isdet = f/ef>=5
    flare_start = min(dt[isdet])-0.01
    flare_end = max(dt[isdet])+0.02
    choose = np.logical_and(dt>=flare_start, dt<=flare_end)
    x = dt[choose]
    y = f[choose]
    ey = ef[choose]
    return x,y,ey


def get_rband_all():
    dat = get_full_opt()
    choose = np.logical_and(
            dat['#instrument']=='TNT/ULTRASPEC', dat['flt']=='r')
    dat = dat[choose]
    f = dat['flux'].values.astype(float)
    ef = dat['unc'].values.astype(float)
    dt = np.array(dat['mjdstart'].values-dat['mjdstart'].values[0]).astype(float)
    return dt, f, ef


def get_rband_noise():
    dat = get_full_opt()
    choose = np.logical_and(
            dat['#instrument']=='TNT/ULTRASPEC', dat['flt']=='r')
    dat = dat[choose]
    f = dat['flux'].values.astype(float)
    ef = dat['unc'].values.astype(float)
    dt = np.array(dat['mjdstart'].values-dat['mjdstart'].values[0]).astype(float)
    choose = dt > 0.08
    x = dt[choose]
    y = f[choose]
    ey = ef[choose]
    return x,y,ey


def plot_ls(ax, x, y, ey, c='k', lab=""):
    """ Plot the Lomb-Scargle Periodogram 

    And the false-alarm power (2.5%)
    """
    ls = LombScargle(x, y, ey)
    frequency, power = ls.autopower(
            method='slow', maximum_frequency=30000) # floating mean
    period_d = 1/frequency
    period_m = period_d*24*60
    ax.plot(period_m, power, c=c, label=lab)
    level = ls.false_alarm_level(0.025, method='bootstrap')
    ax.axhline(y=level, c=c, ls='--', lw=1)
    

if __name__=="__main__":
    figwidth_mm = 183 # Nature standard
    figwidth_in = (figwidth_mm/10)/2.54 # in inches

    # Initialize figure
    fig,axarr = plt.subplots(1,2,figsize=(figwidth_in,figwidth_in*3.5/8))

    # First panel: r-band flare
    ax = axarr[0]
    x,y,ey = get_rband_flare()
    plot_ls(ax, x, y, ey, c=vals.rc, lab='$r$-band flare')
    x,y,ey = get_rband_noise()
    plot_ls(ax, x, y, ey, c='lightgrey', lab='$r$-band noise')
    x,y,ey = get_rband_all()
    plot_ls(ax, x, y, ey, c='k', lab='$r$-band all')
    ax.set_ylabel("Lomb-Scargle Power")
    ax.set_xlabel("Period (Minutes)")
    ax.set_xscale('log')
    ax.set_ylim(0, 0.19)
    ax.set_xlim(0.2, 48) # 48 min = 1/dt
    ax.legend(loc='upper left')

    # First panel: g-band flare
    ax = axarr[1]
    x,y,ey = get_gband_flare()
    plot_ls(ax, x, y, ey, c=vals.gc, lab='$g$-band flare')
    x,y,ey = get_gband_noise()
    plot_ls(ax, x, y, ey, c='lightgrey', lab='$g$-band noise')
    x,y,ey = get_gband_all()
    plot_ls(ax, x, y, ey, c='k', lab='$g$-band all')
    ax.set_ylabel("Lomb-Scargle Power")
    ax.set_xlabel("Period (Minutes)")
    ax.set_xscale('log')
    ax.set_ylim(0, 0.19)
    ax.set_xlim(0.2, 48)
    ax.legend(loc='upper left')

    plt.tight_layout()
    plt.savefig(
            "lomb_scargle_periodogram.eps", dpi=300, 
            pad_inches=0, bbox_inches='tight')
    plt.close()
    #plt.show()
