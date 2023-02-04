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
    choose = dat['Filt']==filt
    x = dat['MJD'][choose].values
    y = dat['Flux'][choose].values
    ey = dat['Unc'][choose].values

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


def plot_keck(ax, dat, t0):
    """ Plot a light-curve panel """
    t0_mjd = float(Time(t0, format='isot').mjd)
    m = 5

    # Get the relevant data
    choose = dat['Filt']=='u'
    x = dat['MJD'][choose].values
    y = dat['Flux'][choose].values
    ey = dat['Unc'][choose].values

    # Plot in minutes
    dt = (x-t0_mjd)*24*60

    # Plot detections
    ax.errorbar(dt, y, ey, c=vals.uc, fmt='*', ms=m*2, lw=0.5, label='$u$')

    # Get the relevant data
    choose = dat['Filt']=='i'
    x = dat['MJD'][choose].values
    y = dat['Flux'][choose].values
    ey = dat['Unc'][choose].values

    # Plot in minutes
    dt = (x-t0_mjd)*24*60

    # Plot detections
    ax.errorbar(dt, y, ey, c=vals.ic, fmt='D', ms=m, lw=0.5, label='$i$')

    # Format x-axis
    t0_str = t0[0:16].replace('T', ' ')
    ax.set_xlabel("Min. since %s" %t0_str, fontsize=8)

    ax.axhline(y=0, c='grey', lw=0.5)
    ax.legend(loc='upper right', ncol=1, fontsize=8)


def plot_magellan(ax):
    """ Plot the Magellan/IMACS flare in uJy """
    tel,mjd,filt,mag,emag,limmag,flare = get_flares()
    choose = tel=='Magellan'
    y = (10**((mag[choose]-8.9)/(-2.5))) * 1E6
    ey = y*emag[choose]*(np.log(10))/2.5
    t0 = Time("2022-12-15T04:30:00", format='isot').mjd
    ax.errorbar((mjd[choose]-t0)*24*60, y, ey,
                fmt='s', c=vals.gc)
    ax.text(0.02, 0.15, 'IMACS $g$-band', transform=ax.transAxes,
            ha='left', va='top', fontsize=8)
    #ax.text(0.02, 0.88, '2022-12-15', transform=ax.transAxes,
    #        ha='left', va='top', fontsize=8)
    ax.set_xlabel("Min. since 2022-12-15 04:30", fontsize=8)


def plot_lt(ax):
    tel,mjd,filt,mag,emag,limmag,flare = get_flares()
    t0 = Time("2022-12-16T20:30:00", format='isot').mjd
    choose = np.logical_and.reduce((tel=='LT/IO:O', flare=='*',
                                    np.abs(mjd-t0)<1))
    y = (10**((mag[choose]-8.9)/(-2.5))) * 1E6
    ey = y*emag[choose]*(np.log(10))/2.5

    gdet = np.logical_and(choose, filt=='g')
    ax.errorbar((mjd[gdet]-t0)*24*60, y, ey,
                fmt='s', c=vals.gc)
    ax.text(0.98, 0.98, 'LT $g$-band', transform=ax.transAxes,
            ha='right', va='top', fontsize=8)
    ax.set_xlabel("Min. since 2022-12-16 20:30", fontsize=8)



def plot_ztf(ax, flarenum=1, window=1):
    """ Plot one of the ZTF flares. Can be flarenum=1 or flarenum=2
    or flarenum=3 (there are two r-band flares and one i-band flare) 

    Window represents the days on either side that you show
    """

    # Get the data
    jd,exp,filt,mag,emag,fujy,efujy = get_ipac()

    # The dates of the two flares
    t0s = np.array(
            ['2022-10-04T09:30:00', '2022-10-05T08:10:04.002', 
             '2022-10-29T04:30:00'])
    t0 = Time(t0s[flarenum-1], format='isot').jd

    # Plot the r-band data from that day 
    choose = np.logical_and.reduce((filt=='r', jd>t0-window, jd<t0+window))
    if sum(choose)>0:
        ax.errorbar((jd[choose]-t0), fujy[choose], efujy[choose], 
                fmt='o', c=vals.rc, label='$r$')

    # Plot the g-band data from that day 
    choose = np.logical_and.reduce((filt=='g', jd>t0-window, jd<t0+window))
    if sum(choose)>0:
        ax.errorbar((jd[choose]-t0), fujy[choose], efujy[choose], 
                    fmt='s', c=vals.gc, label='$g$')

    # Plot the i-band data from that day
    choose = np.logical_and.reduce((filt=='i', jd>t0-window, jd<t0+window))
    if sum(choose)>0:
        ax.errorbar((jd[choose]-t0), fujy[choose], efujy[choose], 
                    fmt='D', c=vals.ic, label='$i$')

    # Plot LT points for that first flare
    if flarenum==1:
        dat = get_keck_lt_ultraspec()
        tel = dat['Tel'].values
        jd = Time(dat['MJD'].values, format='mjd').jd
        choose = np.logical_and.reduce((
            tel=='LT/IO:O', dat['Filt']=='g', jd>t0-window, jd<t0+window))
        x = jd[choose]
        y = dat['Flux'][choose].values
        ey = dat['Unc'][choose].values
        ax.errorbar((x-t0), y, ey, fmt='s', c=vals.gc)

        choose = np.logical_and.reduce((
            tel=='LT/IO:O', dat['Filt']=='r', jd>t0-window, jd<t0+window))
        x = jd[choose]
        y = dat['Flux'][choose].values
        ey = dat['Unc'][choose].values
        ax.errorbar((x-t0), y, ey, fmt='o', c=vals.rc)

    # Formatting the panel
    t0_str = t0s[flarenum-1][0:16].replace('T', ' ')
    return t0_str



if __name__=="__main__":
    # Initialize figure
    fig,axarr = plt.subplots(figsize=(7,7))

    # Plot ZTF: the first r and i flares (1d apart)
    # And also the LT detections of what is presumably the transient
    ax = plt.subplot(4,2,1)
    plot_ztf(ax, window=2.5)
    ax.set_ylabel(r"$f_\nu$ ($\mu$Jy)")
    ax.axhline(y=0, c='grey', lw=0.5)
    axins = ax.inset_axes([0.2, 0.54, 0.2, 0.42])
    t0_str = plot_ztf(axins, window=0.1)
    ax.legend(loc='upper right', fontsize=8, ncol=1, columnspacing=0.2,
          handletextpad=0.1)
    ax.set_xlabel("Days since %s" %t0_str, fontsize=8)
    ax.text(0.01, 0.95, 'ZTF+LT', ha='left', va='top', fontsize=8,
            transform=ax.transAxes)
    axins.set_yticks([])
    axins.set_xticks([])
    axins.set_xlabel("15 min", fontsize=8, labelpad=1)
    ax.indicate_inset_zoom(axins, edgecolor="grey")

    # Plot ZTF: the second r-band flare
    ax = plt.subplot(4,2,2)
    t0_str = plot_ztf(ax, flarenum=3, window=1.5)
    ax.legend(loc='lower center', fontsize=8, ncol=3, columnspacing=0.2,
          handletextpad=0.1)
    ax.set_xlabel("Days since %s" %t0_str, fontsize=8)
    ax.text(0.05, 0.95, 'ZTF', ha='left', va='top', fontsize=8,
            transform=ax.transAxes)
    axins = ax.inset_axes([0.7, 0.54, 0.28, 0.45])
    plot_ztf(axins, flarenum=3, window=0.1)
    axins.set_yticks([])
    axins.set_xticks([])
    axins.set_xlabel("7 min", fontsize=8, labelpad=1)
    ax.indicate_inset_zoom(axins, edgecolor="grey")

    # Plot IMACS panel
    ax = plt.subplot(4,4,5)
    plot_magellan(ax)
    ax.set_ylabel(r"$f_\nu$ ($\mu$Jy)")

    # Plot LT panel
    ax = plt.subplot(4,4,6)
    plot_lt(ax)

    # Plot Keck panel
    ax = plt.subplot(4,4,7)
    dat = get_ultraspec()
    plot_keck(ax, dat, '2022-12-29T10:00:00')
    ax.set_xlim(8, 33)
    ax.text(0.02, 0.02, 'LRIS', transform=ax.transAxes,
            ha='left', va='bottom', fontsize=8)

    # Plot ULTRASPEC r-band panel
    ax = plt.subplot(4,1,3)
    dat = get_ultraspec()
    t0 = Time("2022-12-19T15:00:00", format='isot').mjd
    plot_ultraspec_panel(ax, dat, t0, 'r', 'o', vals.rc)
    ax.set_xlabel("Hours since 2022-12-19 15:00")
    ax.text(0.02, 0.95, 'ULTRASPEC $r$-band', transform=ax.transAxes,
            ha='left', va='top', fontsize=8)
    ax.set_xlim(-0.6, 4.3)
    ax.set_ylim(-6, 33)
    ax.set_ylabel(r"$f_\nu$ ($\mu$Jy)")

    # Zoom-in to r-band flare
    axins = ax.inset_axes([0.5, 0.54, 0.48, 0.45])
    plot_ultraspec_panel(axins, dat, t0, 'r', 'o', vals.rc)
    axins.tick_params(axis='both', labelsize=8, pad=0.5)
    axins.set_ylabel("")
    axins.set_xlim(0.6,1.1)
    axins.set_ylim(-1,32)
    ax.indicate_inset_zoom(axins, edgecolor="grey")

    # Plot ULTRASPEC g-band panel
    ax = plt.subplot(4,1,4)
    t0 = Time("2022-12-20T15:00:00", format='isot').mjd
    plot_ultraspec_panel(ax, dat, t0, 'g', 's', vals.gc)#, plot_binned=True)
    ax.text(0.98, 0.06, 'ULTRASPEC $g$-band', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=8)
    ax.set_xlabel("Time since 2022-12-20 15:00 (hours)")
    ax.set_ylabel(r"$f_\nu$ ($\mu$Jy)")
    ax.set_xlim(0.2, 4.3)

    # Zoom-in to g-band flare
    axins = ax.inset_axes([0.01, 0.54, 0.4, 0.42])
    plot_ultraspec_panel(axins, dat, t0, 'g', 's', vals.gc, plot_binned=True)
    axins.tick_params(axis='both', labelsize=8, pad=0.5)
    axins.set_ylabel("")
    axins.set_xlim(2.5,4.2)
    axins.set_ylim(-1,21)
    ax.indicate_inset_zoom(axins, edgecolor="grey")
    axins.set_xticks([])
    axins.set_yticks([])

    plt.tight_layout()
    plt.show()
    #plt.savefig("flares.png", dpi=300, 
    #            bbox_inches='tight', pad_inches=0.1)
    #plt.close()


