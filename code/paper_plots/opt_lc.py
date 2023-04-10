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
from opt_lc_comparison import *
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


def plot_spec_epochs(ax):
    """ Plot epochs of Keck/LRIS spectroscopy """
    start = 2459846.07834 # from reading the file from Alex's group
    plot_epoch(ax, start, '', '', label=None, c='k', lw=2, ymin=0, ymax=0.05)

    # I think the second one was at 10am ET on Thursday Oct 6,
    # from looking at my email
    start = Time("2022-10-06T15:00:00", format='isot').jd
    plot_epoch(ax, start, '', '', c='k', lw=2, ymin=0, ymax=0.05)


def plot_xray_epochs(ax):
    """ Plot timing of X-ray observations """
    xray = load_swift()
    mjd = xray['!MJD    '].values
    t = Time(mjd, format='mjd').jd

    ch = load_chandra()
    tch = np.array([Time(val, format='isot').jd for val in ch['Start'].values])

    tboth = np.hstack((t, tch))

    for i,tval in enumerate(tboth):
        plot_epoch(
                ax, tval, 'X-ray', align='center', c='k', l=':', lw=1,
                ymin=0.1, ymax=0.15)


def plot_radio_epochs(ax):
    dat = get_radio() 
    xs = np.array([Time(val, format='isot').jd for val in dat['Date'].values])
    for i,x in enumerate(xs):
        plot_epoch(ax, x, 'Radio', align='right', ymin=0.05, ymax=0.1, lw=1)


def plot_flare_epochs(ax, dat):
    jd = Time(dat['mjdstart'].values, format='mjd').jd
    isflare = dat['isflare'].values
    filts = dat['flt'].values
    cs = [vals.rc, vals.gc, vals.ic, vals.uc, vals.wc]

    for j,filt in enumerate(np.array(['r', 'g', 'i', 'u', 'w'])):
        choose = np.logical_and(isflare, filts==filt)
        flare_epochs = jd[choose]
        for i,x in enumerate(flare_epochs):
            plot_epoch(ax, x, 'Opt. Flare', align='right', 
                       ymin=0.95, ymax=1, c=cs[j], lw=1)


if __name__=="__main__":
    """ Plot the FULL light curve, with all detections """
    # Get data
    dat = get_full_opt()

    # Get the approximate time of explosion
    t0_str = Time(vals.t0, format='jd').isot.replace("T", " ").split('.')[0]

    # Get the time where you split the axes
    tsplit = 20

    # Initialize
    fig,axarr = plt.subplots(
            1,2,figsize=(7,3),sharey=True,gridspec_kw={'width_ratios': [1.5,3]})

    # AT2022tsd in both panels
    for ax in axarr:
        plot_22tsd(ax, show='apparent', offset=0.11)

    # Left panel: comparisons
    ax = axarr[0]
    shift = np.abs(Planck18.distmod(z=0.0141).value-vals.dm)
    plot_18cow(ax, show='apparent', offset=shift)
    shift = np.abs(Planck18.distmod(z=0.2442).value-vals.dm)
    plot_20xnd(ax, show='apparent', offset=shift)
    shift = np.abs(Planck18.distmod(z=0.0085).value-vals.dm)
    plot_98bw(ax, show='apparent', offset=shift)
    shift = np.abs(Planck18.distmod(z=0.677).value-vals.dm)
    plot_sn2011kl(ax, show='apparent', offset=shift)

    # Left panel: spec epochs
    plot_spec_epochs(ax)

    # Left panel: formatting
    ax.set_yticks([19,20,21,22,23])
    ax.set_yticklabels([19,20,21,22,23])
    ax.set_ylabel(r"$m_\mathrm{opt}$ (AB)", fontsize=10,
            fontname='sans-serif')
    ax.set_xlim(-2,tsplit)

    # Right panel: comparisons and epochs
    ax = axarr[1]
    plot_spec_epochs(ax)
    ax.text(15.4, 22.95, 'Opt. spec.', fontsize=8, ha='right', c='k')
    plot_xray_epochs(ax)
    ax.text(24, 22.75, 'X-ray', fontsize=8, ha='right', c='k')
    plot_radio_epochs(ax)
    ax.text(26.2, 22.55, 'Radio', fontsize=8, ha='right', c='k')
    plot_flare_epochs(ax, dat)
    ax.text(26.2, 19.2, 'Opt. flares', fontsize=8, ha='right', c='k')
    shift = np.abs(Planck18.distmod(z=0.0085).value-vals.dm)
    plot_98bw(ax, show='apparent', offset=shift)
    shift = np.abs(Planck18.distmod(z=0.677).value-vals.dm)
    plot_sn2011kl(ax, show='apparent', offset=shift)

    # Formatting of the right axis
    ax.set_xlim(tsplit, 210)
    ax.set_xscale('log')
    ax.set_ylim(23, 19)
    ax.set_xticks([30,40,50,70,100,200])
    ax.set_xticklabels([30,40,50,70,100,200])

    # Make a second y-axis
    ax2 = ax.twinx()
    ax2.set_ylabel(
            "$M_\mathrm{opt} (AB)$", fontsize=10, rotation=270, labelpad=15.0)
    y_f = lambda y_i: y_i-vals.dm+2.5*np.log10(1+vals.z)
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.tick_params(axis='y', labelsize=10)
    ax2.set_yticks([-21,-20,-19,-18])
    ax2.set_yticklabels([-21,-20,-19,-18])
    ax2.plot([],[])

    # For the legend
    ax.plot([100,100],[150,150],ls='-.',c='grey',label='SN2011kl')
    ax.plot([100,100],[150,150],ls='--',c=vals.gc,label='AT2018cow $g$',lw=0.5)
    ax.plot([100,100],[150,150],ls='--',c=vals.rc,label='AT2018cow $r$',lw=0.5)
    ax.plot([100,100],[150,150],ls='-',c=vals.gc,label='AT2020xnd $g$',lw=0.5)
    ax.plot([100,100],[150,150],ls='-',c=vals.rc,label='AT2020xnd $r$',lw=0.5)
    ax.plot([100,100],[150,150],ls=':',c=vals.rc,label='SN1998bw')
    ax.scatter(0,0,marker='s',c=vals.gc, label='AT2022tsd $g$')
    ax.scatter(0,0,marker='o',c=vals.rc, label='AT2022tsd $r$')
    ax.scatter(0,0,marker='D',c=vals.ic, label='AT2022tsd $i$')
    fig.legend(fontsize=8, bbox_to_anchor=(0.5,1.02), loc='upper center',
               ncol=5, handletextpad=0.1)

    # Formatting of the whole fig
    fig.text(0.5, 0, r'$\Delta t_\mathrm{obs}$ (days)', ha='center')
    axarr[0].spines['right'].set_visible(False)
    axarr[1].spines['left'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    axarr[0].yaxis.tick_left()
    axarr[0].tick_params(labelright='off')
    axarr[1].yaxis.tick_right()
    axarr[1].set_yticks([])
    ax2.yaxis.tick_right()

    fig.subplots_adjust(wspace=0)
    #plt.show()
    plt.savefig("opt_lc.png", dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close()
