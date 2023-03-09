""" Plot the full optical light curve, with zoom-ins on the discovery period
and on the flares. """

import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.patches as patches
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["pdf.fonttype"] = 42
from astropy.time import Time
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
from get_opt import *
from get_radio_at2022tsd import *
from get_xray import *
import vals


def plot_lim(ax, x, y, band, leg=False):
    """ Plot 5-sigma upper limits """
    col = vals.gc
    lw = 0.3
    s = 10
    lab = None
    if band=='r':
        col = vals.rc
    elif band=='w':
        col = vals.wc
    if leg:
        lab = '$%s$ limit' %band
    ax.scatter(x, y, 
        edgecolor=col, facecolor='white', marker='v', lw=lw, s=s, label=lab)


def plot_det(ax, x, y, ey, band, lines=False, leg=False):
    """ Plot detections """
    col = vals.gc
    m = 's'
    s = 4
    lab = None # Default is no label
    if band=='r':
        col = vals.rc
        m = 'o'
    elif band=='i':
        col = vals.ic
        m = 'D'
    elif band=='u':
        col = vals.uc
        m = '>'
    elif band=='w':
        col = vals.wc
        m = '<'
    if lines:
        m = '%s-'%m
    if leg:
        lab='$%s$' %band
    ax.errorbar(x, y, ey, c=col, fmt=m, label=lab, ms=s, lw=0.5)


def plot_full_lc(ax, dat):
    """ Add all the optical photometry """
    jd = Time(dat['mjdstart'].values, format='mjd').jd
    dt = jd-vals.t0

    isflare = dat['isflare']
    istransient = dat['istransient']
    mag = dat['mag_extcorr']
    limmag = dat['maglim_extcorr']
    emag = dat['emag']
    filt = dat['flt']

    isdet = np.logical_or(isflare, istransient)

    # Plot the g-band detections
    choose = np.logical_and(isdet, filt=='g')
    plot_det(ax,dt[choose],mag[choose],emag[choose],'g', leg=True)

    # Plot the g-band limits
    choose = np.logical_and(~isdet, filt=='g')
    plot_lim(ax, dt[choose], limmag[choose], 'g', leg=True)

    # Now onto r-band
    # HCT is technically R band
    dat.loc[dat['flt']=='R', 'flt'] = 'r'

    # Plot the r-band flares
    choose = np.logical_and(isdet, filt=='r')
    plot_det(ax,dt[choose],mag[choose],emag[choose],'r', leg=True)

    # Plot the r-band limits
    choose = np.logical_and(filt=='r', ~isdet)
    plot_lim(ax, dt[choose], limmag[choose], 'r', leg=True)

    # Plot the i-band flares ( will include PS 1)
    choose = np.logical_and(isdet, filt=='i')
    plot_det(ax,dt[choose],mag[choose],emag[choose],'i',
             leg=True)

    # Plot the i-band limits
    choose = np.logical_and(filt=='i', ~isdet)
    plot_lim(ax, dt[choose], limmag[choose], 'i', leg=True)

    # Finally, w-band
    choose = np.logical_and(isdet, filt=='w')
    plot_det(ax,dt[choose],mag[choose],emag[choose],'w',
             leg=True)
    choose = np.logical_and(filt=='w', ~isdet)
    plot_lim(ax, dt[choose], limmag[choose], 'w', leg=True)


def plot_flares_zoom(ax):
    """ Plot a zoom-in of the flares. Minutes timescale """
    # Non-ZTF flares
    tel,mjd,filt,mag,emag,limmag,flare = get_flares()
    jd = Time(mjd, format='mjd').jd

    # The Magellan flare is in g-band
    choose = np.logical_and(flare=='*', tel=='Magellan')
    x = (jd[choose]-vals.t0)*24*60
    y = mag[choose]-vals.ext['g']
    ey = emag[choose]

    # Convert to luminosity
    # Frequency: Sloan g' for IMACS
    freq = 3E18/(4723.59) 
    fjy = 3631*10**(y/(-2.5))
    Lnu = fjy * 1E-23 * 4 * np.pi * (vals.dL_cm)**2
    vLv = Lnu * freq 
    ey = np.zeros(len(vLv))

    plot_det(ax, x-x[0], vLv, ey, 'g')


def plot_epoch(ax, xval, txt, align='center', label=None, ymin=0, ymax=0.05, c='grey', l='-', lw=2):
    """ Plot an epoch: a vertical line at x, labeled with txt 
    xval should be in JD """
    ax.axvline(
            x=xval-vals.t0, ymin=ymin, ymax=ymax, 
            c=c, ls=l, lw=lw, label=label)
    #ax.text(
    #        xval-vals.t0, 23.9, txt, ha=align, va='bottom',
    #        color='grey', rotation=270, fontsize=8)


def plot_spec_epochs(ax):
    """ Plot epochs of Keck/LRIS spectroscopy """
    start = 2459846.07834 # from reading the file from Alex's group
    plot_epoch(ax, start, '', '', label=None, c='k', lw=2)

    # I think the second one was at 10am ET on Thursday Oct 6,
    # from looking at my email
    start = Time("2022-10-06T15:00:00", format='isot').jd
    plot_epoch(ax, start, '', '', label='Spec.', c='k', lw=2)


def plot_xray_epochs(ax):
    """ Plot timing of X-ray observations """
    xray = load_swift()
    mjd = xray['!MJD    '].values
    t = Time(mjd, format='mjd').jd

    ch = load_chandra()
    tch = np.array([Time(val, format='isot').jd for val in ch['Start'].values])

    tboth = np.hstack((t, tch))

    for i,tval in enumerate(tboth):
        if i==0:
            plot_epoch(
                    ax, tval, 'X-ray', align='center', label='X-ray', 
                    c='k', l=':', lw=1, ymin=0.1, ymax=0.15)
        else:
            plot_epoch(
                    ax, tval, 'X-ray', align='center', c='k', l=':', lw=1,
                    ymin=0.1, ymax=0.15)


def plot_radio_epochs(ax):
    dat = get_radio() 
    xs = np.array([Time(val, format='isot').jd for val in dat['Date'].values])
    for i,x in enumerate(xs):
        print(x)
        if i==0:
            plot_epoch(ax, x, 'Radio', align='right', label='Radio', 
                       ymin=0.05, ymax=0.1, lw=1)
        else:
            plot_epoch(ax, x, 'Radio', align='right', ymin=0.05, ymax=0.1, lw=1)


if __name__=="__main__":
    # Get the full LC 
    dat = get_full_opt()

    # Get the approximate time of explosion
    t0_str = Time(vals.t0, format='jd').isot.replace("T", " ").split('.')[0]

    # Initialize
    fig,ax = plt.subplots(1,1,figsize=(6,3.5))

    # Plot LC 
    plot_full_lc(ax, dat)

    # Plot epochs of various things (spec, x-ray, radio)
    plot_spec_epochs(ax)
    plot_xray_epochs(ax)
    plot_radio_epochs(ax)

    # Make a legend
    #ax.legend(loc='lower right', ncol=1, fontsize=8)

    # Zoom in
    #axins = ax.inset_axes([0.35, 0.15, 0.3, 0.35])
    #plot_flares_zoom(axins)
    #axins.set_xlabel("Minutes (obs. frame)", fontsize=8, labelpad=1)
    #axins.set_ylabel(r"$\nu L_\nu$ (erg s$^{-1}$)", fontsize=8, labelpad=1)
    #axins.tick_params(axis='both', labelsize=8, pad=0.5)
    #axins.set_yscale('log')
    #axins.set_ylim(3E42, 1E44)
    # rect = patches.Rectangle((98.2, 19.4), 3.6, 3.6, linewidth=0.5, 
    #                          edgecolor='grey', facecolor='none')
    # ax.add_patch(rect)
    # ax.plot([92.3,98.2],[21.2,19.4], c='grey', lw=0.5)
    # ax.plot([92.3,98.2],[22.8,19.4+3.6], c='grey', lw=0.5)
    # # Make a second x-axis
    # axins2 = axins.twiny()
    # axins2.set_xlabel(r"Minutes (rest frame)", fontsize=8, labelpad=1)
    # x_f = lambda x_i: x_i/(1+float(vals.z))
    # xmin, xmax = axins.get_xlim()
    # axins2.set_xlim((x_f(xmin), x_f(xmax)))
    # axins2.tick_params(axis='both', labelsize=8, pad=0.5)

    # Formatting
    ax.set_xlim(-5, 145)
    ax.set_ylim(24.5, 18.7)

    # Make a second x-axis
    # ax3 = ax.twiny()
    # ax3.set_xlabel(r"Days since %s (rest frame)" %t0_str)
    # x_f = lambda x_i: x_i/(1+float(vals.z))
    # xmin, xmax = ax.get_xlim()
    # ax3.set_xlim((x_f(xmin), x_f(xmax)))
    # ax3.tick_params(axis='x', labelsize=10)
    # ax3.plot([],[])

    ax.set_xlabel(
            r"Days since %s (observer frame)" %t0_str,fontsize=10,
            fontname='sans-serif')

    ax.set_ylabel(r"Apparent Magnitude", fontsize=10,
            fontname='sans-serif')

    # Make a second y-axis
    ax2 = ax.twinx()
    ax2.set_ylabel("Absolute Magnitude", fontsize=10, rotation=270, labelpad=15.0)
    y_f = lambda y_i: y_i-vals.dm+2.5*np.log10(1+vals.z)
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.tick_params(axis='y', labelsize=10)
    ax2.plot([],[])

    ax.legend(loc='upper center', fontsize=7, handletextpad=0.1, 
              bbox_to_anchor=(0.5, 1.15), 
              ncol=6, fancybox=True)

    #plt.tight_layout()
    #plt.show()
    plt.savefig("opt_lc.png", dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close()
