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


def plot_ultraspec_panel(ax, dat, t0, filt, m, col, plot_binned=False, y2=False):
    """ Plot a light-curve panel """

    # Get ULTRASPEC data
    choose = dat['flt']==filt
    x = dat['mjdstart'][choose].values
    y = dat['flux_extcorr'][choose].values.astype(float)
    ey = dat['unc_extcorr'][choose].values.astype(float)

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
        x = dt[bsize:len(dt)-bsize]
        y = np.array(binned)
        order = np.argsort(x)
        ax.plot(x[order], y[order], c='k', zorder=10, lw=0.5)
        
    if y2:
        # Make an axis on the right-hand side with vLv
        ax2 = ax.twinx()
        leff = vals.sdss_pivot['g'] # use g-band for all
        nueff = 3E18/leff
        y_f = lambda y_i: nueff * y_i * 1E-6 * 1E-23 * \
                4 * np.pi * vals.dL_cm**2 / 1E43
        ymin, ymax = ax.get_ylim()
        ax2.set_ylim((y_f(ymin), y_f(ymax)))
        ax2.plot([],[])
        #ax2.set_yscale('log')
        ax2.tick_params(axis='both', labelsize=9)
        ax2.set_ylabel(
                r"$\nu L_\nu$ ($10^{43}$ erg s$^{-1}$)", fontsize=9, 
                rotation=270, labelpad=15.0)


def plot_flare(ax, tab, mjd, window=1, unit='Days'):
    """ Plot the flare from a single night 
    The window is the interval used to select data around the flare. 

    Make sure to correct for Galactic extinction
    And then plot vLv on the other y-axis
    """

    # Data from the night of choice
    choose = np.logical_and(tab['flare_nights']==mjd, tab['isflare'])
    indmax = np.argmax(tab['flux'][choose])
    t0 = tab['mjdstart'][choose].values[indmax]

    # Extract relevant data
    filt = tab['flt'].values
    t = tab['mjdstart'].values
    fujy = tab['flux_extcorr'].values
    efujy = tab['unc_extcorr'].values

    # Telescope string
    tel_list = []

    # Plot
    cols = [vals.rc, vals.gc, vals.ic, vals.uc, vals.wc]
    syms = ['o', 's', 'D', '*', '<']
    ms_scale = [1,1,1,2,1.5]
    done = False
    for i,b in enumerate(np.array(['r', 'g', 'i', 'u', 'w'])):
        choose = np.logical_and.reduce((filt==b, t>t0-window, t<t0+window))
        if sum(choose)>0:
            for val in np.unique(tab['#instrument'].values[choose]):
                tel_list.append(val.split('/')[0])
            fac = 1
            if unit=='Hours':
                fac = 24
            if unit=='Minutes':
                fac = 24*60
            ax.errorbar((t[choose]-t0)*fac, fujy[choose], efujy[choose], 
                    fmt=syms[i], c=cols[i], label='$%s$' %b, ms=5*ms_scale[i])

    t0_str = str(Time(t0, format='mjd').isot).replace('T', ' ')[0:-7]
    tel_list = np.unique(np.array(tel_list))
    tel_str = '/'.join(tel_list)

    # Telecope string
    return t0_str,tel_str


if __name__=="__main__":
    # Initialize figure
    fig,axarr = plt.subplots(figsize=(7,8))

    # Get the optical photometry
    tab = get_full_opt()

    # Identify nights with flares (since we'll plot each night individually)
    flare_nights = np.unique(tab['mjdstart'][tab['isflare']==True].astype(int))
    tab['flare_nights'] = tab['mjdstart'].astype(int)

    # There are 12 flare nights. 

    # Some custom settings
    windows = (np.ones(len(flare_nights))/24/60)*30
    windows[0] = 2.5 # first ZTF night needs longer
    units = ['Minutes']*len(flare_nights)
    units[0] = 'Days' # for that first ZTF obs

    for i,night in enumerate(flare_nights):
        if i not in [1,7,8]: # skip three: ZTF i-band, 2xTNT
            if i < 1:
                ind = i + 1
            elif i < 9:
                ind = i
            else:
                ind = i-2
            ax = plt.subplot(5,3,ind)
            t0_str,tel_str = plot_flare(ax, tab, night, window=windows[i],
                                        unit=units[i])
            if i == 0:
                ax.text(0.01, 0.95, tel_str, ha='left', va='top', fontsize=8,
                        transform=ax.transAxes)
            else:
                ax.text(0.99, 0.95, tel_str, ha='right', va='top', fontsize=8,
                        transform=ax.transAxes)
            unit_str = units[i]
            if unit_str=='Minutes':
                unit_str='Min.'
            ax.set_xlabel("%s since %s" %(unit_str,t0_str[5:]), fontsize=8)

            if i==0:
                axins = ax.inset_axes([0.15, 0.45, 0.25, 0.3])
                t0_str,tel_str = plot_flare(axins, tab, night, window=0.2)
                axins.set_yticks([])
                axins.set_xticks([])
                axins.set_xlabel("15 min", fontsize=8, labelpad=1)
                ax.indicate_inset_zoom(axins, edgecolor="grey")
                ax.set_ylabel(r"$f_\nu$ ($\mu$Jy)")
            if i in [0,5,6,7,8,9]:
                ax.axhline(y=0, c='grey', lw=0.5)
            if np.logical_or(tel_str=='NOT', tel_str=='LT+NOT+ZTF'):
                ax.set_ylabel(r"$f_\nu$ ($\mu$Jy)")
            ax.tick_params(axis='both', labelsize=9, pad=0.5)

            # Make a second axis
            scale_y2 = 42
            if ind < 7:
                scale_y2 = 43
            ax2 = ax.twinx()
            leff = vals.sdss_pivot['g'] # use g-band for all
            nueff = 3E18/leff
            y_f = lambda y_i: nueff * y_i * 1E-6 * 1E-23 * \
                    4 * np.pi * vals.dL_cm**2 / 10**scale_y2
            ymin, ymax = ax.get_ylim()
            ax2.set_ylim((y_f(ymin), y_f(ymax)))
            ax2.plot([],[])
            #ax2.set_yscale('log')
            ax2.tick_params(axis='both', labelsize=9)
            if np.logical_or(ind==4, ind==7):
                ax2.set_ylabel(
                        r"$\nu L_\nu$ ($10^{%s}$ erg s$^{-1}$)" %scale_y2, fontsize=9,
                        rotation=270, labelpad=15.0)

            if ind==2:
                ax.set_yticks([5,7,10])
                ax.set_yticklabels([5,7,10])
                ax2.set_yticks([0.7,1.1])
                ax2.set_yticklabels([0.7,1.1])

    # # Plot ULTRASPEC r-band panel
    ax = plt.subplot(5,1,4)
    t0 = Time("2022-12-19T15:00:00", format='isot').mjd
    plot_ultraspec_panel(ax, tab, t0, 'r', 'o', vals.rc, y2=True)
    ax.set_xlabel("Hours since 12-19 15:00", fontsize=9)
    ax.text(0.02, 0.95, 'ULTRASPEC $r$-band', transform=ax.transAxes,
            ha='left', va='top', fontsize=8)
    ax.set_xlim(-0.6, 4.3)
    ax.set_ylim(-6, 33)
    ax.set_ylabel(r"$f_\nu$ ($\mu$Jy)")
    ax.tick_params(axis='both', labelsize=9, pad=0.5)

    # Zoom-in to r-band flare
    axins = ax.inset_axes([0.5, 0.54, 0.48, 0.45])
    plot_ultraspec_panel(axins, tab, t0, 'r', 'o', vals.rc)
    axins.tick_params(axis='both', labelsize=8, pad=0.5)
    axins.set_ylabel("")
    axins.set_xlim(0.6,1.1)
    axins.set_ylim(-1,32)
    ax.indicate_inset_zoom(axins, edgecolor="grey")

    # Plot ULTRASPEC g-band panel
    ax = plt.subplot(5,1,5)
    t0 = Time("2022-12-20T15:00:00", format='isot').mjd
    plot_ultraspec_panel(ax, tab[tab['#instrument']=='TNT/ULTRASPEC'], 
                         t0, 'g', 's', vals.gc, y2=True)
    ax.text(0.02, 0.95, 'ULTRASPEC $g$-band', transform=ax.transAxes,
            ha='left', va='top', fontsize=8)
    ax.set_xlabel("Hours since 12-20 15:00", fontsize=9)
    ax.set_ylabel(r"$f_\nu$ ($\mu$Jy)")
    ax.set_xlim(0.2, 4.3)
    ax.set_ylim(-5,22)

    # Zoom-in to g-band flare
    axins = ax.inset_axes([0.01, 0.35, 0.4, 0.42])
    plot_ultraspec_panel(axins, tab[tab['#instrument']=='TNT/ULTRASPEC'], 
                         t0, 'g', 's', vals.gc, plot_binned=True)
    axins.tick_params(axis='both', labelsize=8, pad=0.5)
    axins.set_ylabel("")
    axins.set_xlim(2.5,4.2)
    axins.set_ylim(-1,21)
    ax.indicate_inset_zoom(axins, edgecolor="grey")
    axins.set_xticks([])
    axins.set_yticks([])

    # Add an axis just for the legend?
    #ax = plt.subplot(4,4,8)
    #ax.scatter(0,0,c=vals.uc,marker='*',label='$u$',s=60)
    #ax.scatter(0,0,c=vals.gc,marker='s',label='$g$')
    #ax.scatter(0,0,c=vals.rc,marker='o',label='$r$')
    #ax.scatter(0,0,c=vals.ic,marker='D',label='$i$')
    #ax.set_xlim(5,10)
    #ax.axis('off')
    #ax.legend()

    plt.subplots_adjust(wspace=0.4, hspace=0.5)
    plt.show()
    #plt.savefig("flares.png", dpi=300, 
    #            bbox_inches='tight', pad_inches=0.1)
    #plt.close()


