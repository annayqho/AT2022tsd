""" Plot the ULTRASPEC flares """

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import matplotlib.patches as patches
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["pdf.fonttype"] = 42
from astropy.time import Time
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
from get_opt import *
import vals


def plot_ultraspec_panel(ax, dat, t0, filt, m, col, plot_binned=False):
    """ Plot a light-curve panel """

    # Get ULTRASPEC data
    choose = dat['flt']==filt
    x = dat['mjdstart'][choose].values
    y = dat['flux'][choose].values
    ey = dat['unc'][choose].values

    # Marker formatting
    s=2
    l=0.2

    # Definition of a detection: 5-sigma
    # Otherwise you get what seem to be spurious points
    det = y/ey >= 5

    # Plot in hours
    dt = (x-t0)*24

    # Plot detections
    ax.errorbar(dt[det], y[det], ey[det], c=col, fmt=m, ms=s,lw=l)

    # Plot the non-detections
    ax.errorbar(dt[~det], y[~det], ey[~det],lw=l,
                mec=col, fmt=m, ms=s, alpha=0.3, mfc='white', c=col)

    # Plot a y=0 line
    ax.axhline(y=0, c='grey', lw=0.5)

    # Bin the light curve...group by five points
    if plot_binned:
        bsize = 3
        binned = []
        for i in np.arange(bsize,len(dt)-bsize):
            binned.append(np.average(y[i-bsize:i+bsize], 
                          weights=1/(ey[i-bsize:i+bsize])**2))
        ax.plot(dt[bsize:len(dt)-bsize], np.array(binned), c='k', zorder=10,
                lw=0.5)


def plot_flare(ax, tab, mjd, window=1, unit='Days'):
    """ Plot the flare from a single night 
    The window is the interval used to select data around the flare. """
    # Data from the night of choice
    choose = np.logical_and(tab['flare_nights']==mjd, tab['flare'])
    indmax = np.argmax(tab['flux'][choose])
    t0 = tab['mjdstart'][choose].values[indmax]

    # Extract relevant data
    filt = tab['flt'].values
    t = tab['mjdstart'].values
    fujy = tab['flux'].values
    efujy = tab['unc'].values

    # Telescope string
    tel_list = []

    # Plot
    cols = [vals.rc, vals.gc, vals.ic, vals.uc]
    syms = ['o', 's', 'D', '*']
    for i,b in enumerate(np.array(['r', 'g', 'i', 'u'])):
        choose = np.logical_and.reduce((filt==b, t>t0-window, t<t0+window))
        if sum(choose)>0:
            for val in np.unique(tab['#inst'].values[choose]):
                tel_list.append(val.split('/')[0])
            fac = 1
            if unit=='Hours':
                fac = 24
            if unit=='Minutes':
                fac = 24*60
            ax.errorbar((t[choose]-t0)*fac, fujy[choose], efujy[choose], 
                    fmt=syms[i], c=cols[i], label='$%s$' %b)

    t0_str = str(Time(t0, format='mjd').isot).replace('T', ' ')[0:-7]
    tel_list = np.unique(np.array(tel_list))
    tel_str = '+'.join(tel_list)

    # Telecope string
    return t0_str,tel_str


if __name__=="__main__":
    # Initialize figure
    fig,axarr = plt.subplots(figsize=(7,7))

    # Get the optical photometry
    tab = get_full_opt()

    # Identify flares as 5-sigma detections after a chosen time
    tab['flare'] = np.logical_and(tab['sig']>5, tab['mjdstart']>59856)

    # Identify nights with flares (since we'll plot each night individually)
    flare_nights = np.unique(tab['mjdstart'][tab['flare']].astype(int))
    tab['flare_nights'] = tab['mjdstart'].astype(int)

    # There are 10 flare nights. 

    # Some custom settings
    windows = (np.ones(len(flare_nights))/24/60)*30
    windows[0] = 2.5 # first ZTF night needs longer
    units = ['Minutes']*len(flare_nights)
    units[0] = 'Days' # for that first ZTF obs

    for i,night in enumerate(flare_nights):
        if i != 1: # skip the ZTF i-band point, already in first panel
            if i < 1:
                ind = i + 1
            else:
                ind = i
            ax = plt.subplot(4,3,ind)
            t0_str,tel_str = plot_flare(ax, tab, night, window=windows[i],
                                        unit=units[i])
            ax.text(0.01, 0.95, tel_str, ha='left', va='top', fontsize=8,
                    transform=ax.transAxes)
            ax.set_ylabel(r"$f_\nu$ ($\mu$Jy)")
            ax.axhline(y=0, c='grey', lw=0.5)
            ax.set_xlabel("%s since %s" %(units[i],t0_str), fontsize=8)

            if i==0:
                axins = ax.inset_axes([0.2, 0.4, 0.2, 0.42])
                t0_str,tel_str = plot_flare(axins, tab, night, window=0.2)
                axins.set_yticks([])
                axins.set_xticks([])
                axins.set_xlabel("15 min", fontsize=8, labelpad=1)
                ax.indicate_inset_zoom(axins, edgecolor="grey")

    # # Plot ULTRASPEC r-band panel
    # ax = plt.subplot(4,1,3)
    # dat = get_non_ztf()
    # t0 = Time("2022-12-19T15:00:00", format='isot').mjd
    # plot_ultraspec_panel(ax, dat, t0, 'r', 'o', vals.rc)
    # ax.set_xlabel("Hours since 2022-12-19 15:00")
    # ax.text(0.02, 0.95, 'ULTRASPEC $r$-band', transform=ax.transAxes,
    #         ha='left', va='top', fontsize=8)
    # ax.set_xlim(-0.6, 4.3)
    # ax.set_ylim(-6, 33)
    # ax.set_ylabel(r"$f_\nu$ ($\mu$Jy)")

    # # Zoom-in to r-band flare
    # axins = ax.inset_axes([0.5, 0.54, 0.48, 0.45])
    # plot_ultraspec_panel(axins, dat, t0, 'r', 'o', vals.rc)
    # axins.tick_params(axis='both', labelsize=8, pad=0.5)
    # axins.set_ylabel("")
    # axins.set_xlim(0.6,1.1)
    # axins.set_ylim(-1,32)
    # ax.indicate_inset_zoom(axins, edgecolor="grey")

    # # Plot ULTRASPEC g-band panel
    # ax = plt.subplot(4,1,4)
    # t0 = Time("2022-12-20T15:00:00", format='isot').mjd
    # plot_ultraspec_panel(ax, dat, t0, 'g', 's', vals.gc)#, plot_binned=True)
    # ax.text(0.98, 0.06, 'ULTRASPEC $g$-band', transform=ax.transAxes,
    #         ha='right', va='bottom', fontsize=8)
    # ax.set_xlabel("Time since 2022-12-20 15:00 (hours)")
    # ax.set_ylabel(r"$f_\nu$ ($\mu$Jy)")
    # ax.set_xlim(0.2, 4.3)

    # # Zoom-in to g-band flare
    # axins = ax.inset_axes([0.01, 0.54, 0.4, 0.42])
    # plot_ultraspec_panel(axins, dat, t0, 'g', 's', vals.gc, plot_binned=True)
    # axins.tick_params(axis='both', labelsize=8, pad=0.5)
    # axins.set_ylabel("")
    # axins.set_xlim(2.5,4.2)
    # axins.set_ylim(-1,21)
    # ax.indicate_inset_zoom(axins, edgecolor="grey")
    # axins.set_xticks([])
    # axins.set_yticks([])

    plt.tight_layout()
    plt.show()
    #plt.savefig("flares.png", dpi=300, 
    #            bbox_inches='tight', pad_inches=0.1)
    #plt.close()


