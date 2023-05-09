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
    """ Plot the flare epochs. Use a vertical line stretching from
    the minimum value to the maximum value. """
    jd = Time(dat['mjdstart'].values, format='mjd').jd
    isflare = dat['isflare'].values
    filts = dat['flt'].values
    cs = [vals.rc, vals.gc, vals.ic, vals.uc, vals.wc, vals.wc]

    for j,filt in enumerate(np.array(['r', 'g', 'i', 'u', 'w', 'clear'])):
        choose = np.logical_and(isflare, filts==filt)
        flare_epochs = jd[choose]
        flare_epochs_int = flare_epochs.astype(int)
        flare_epochs_int_unique = np.unique(flare_epochs_int)
        # for each epoch,
        for i,night in enumerate(flare_epochs_int_unique):
            choose_data = flare_epochs_int==night
            mags = dat['mag_extcorr'][choose].values[choose_data]
            emags = dat['emag'][choose].values[choose_data]
            dt_night = night-vals.t0
            if len(mags)>1:
                ax.plot([dt_night,dt_night], [min(mags),max(mags)], 
                        c=cs[j], lw=2)
            else:
                ax.plot([dt_night,dt_night], 
                        [mags[0]-emags[0],mags[0]+emags[0]], c=cs[j], lw=2)


def plot_nonflare_epochs(ax, dat):
    """ Plot the non-flare epochs. Show an upper limit
    corresponding to the median limit or something like that. """
    mjd = dat['mjdstart'].values
    jd = Time(mjd, format='mjd').jd
    isflare = dat['isflare'].values
    istransient = dat['istransient'].values
    nonflare = np.logical_and.reduce((~istransient, ~isflare))
    filts = dat['flt'].values
    cs = [vals.rc, vals.gc, vals.ic, vals.uc, vals.wc]

    for j,filt in enumerate(np.array(['r', 'g', 'i', 'u', 'w'])):
        # Identify nights with NO flares
        has_flare = np.logical_and(isflare, filts==filt)
        night_has_flare = np.unique(jd[has_flare].astype(int))
        choose = np.logical_and(~isflare, filts==filt)
        flare_epochs = jd[choose]
        flare_epochs_int = flare_epochs.astype(int)
        flare_epochs_int_unique = np.unique(flare_epochs_int)
        # for each epoch,
        for i,night in enumerate(flare_epochs_int_unique):
            if night not in night_has_flare:
                choose_data = flare_epochs_int==night
                mags = dat['maglim_extcorr'][choose].values[choose_data]
                dt_night = night-vals.t0
                show_maglim = mags[0]
                if len(mags)>1:
                    show_maglim = np.median(mags)
                ax.scatter(dt_night, show_maglim, edgecolor=cs[j], marker='v',
                           facecolor='white', zorder=0, lw=0.5, s=5)
                          


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
        plot_22tsd(ax, show='apparent')
        plot_nonflare_epochs(ax, dat)

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
    shift = np.abs(Planck18.distmod(z=0.1353).value-vals.dm)
    plot_at2020mrf(ax, show='apparent', offset=shift-1.5)

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
    ax.text(16, 23.45, 'Opt. spec.', fontsize=8, ha='right', c='k')
    plot_xray_epochs(ax)
    ax.text(24.3, 23.25, 'X-ray', fontsize=8, ha='right', c='grey')
    plot_radio_epochs(ax)
    ax.text(26.7, 23, 'Radio', fontsize=8, ha='right', c='k')
    plot_flare_epochs(ax, dat)
    ax.text(26.6, 19.1, '$i$-band flare', fontsize=8, ha='right', c=vals.ic)
    ax.text(25.8, 20.1, '$r$-band flare', fontsize=8, ha='right', c=vals.rc)
    ax.text(68, 21.1, '$w$-band/clear flare', fontsize=8, ha='right', c=vals.wc)
    ax.text(97, 20, '$g$-band flare', fontsize=8, ha='right', c=vals.gc)
    shift = np.abs(Planck18.distmod(z=0.0085).value-vals.dm)
    plot_98bw(ax, show='apparent', offset=shift)
    shift = np.abs(Planck18.distmod(z=0.677).value-vals.dm)
    plot_sn2011kl(ax, show='apparent', offset=shift)
    shift = np.abs(Planck18.distmod(z=0.1353).value-vals.dm)
    plot_at2020mrf(ax, show='apparent', offset=shift)

    # Formatting of the right axis
    ax.set_xlim(tsplit, 210)
    ax.set_xscale('log')
    ax.set_ylim(23.5, 18.8)
    ax.set_xticks([30,40,50,70,100,200])
    ax.set_xticklabels([30,40,50,70,100,200])

    # Make a second y-axis
    ax2 = ax.twinx()
    ax2.set_ylabel(
            "$M_\mathrm{opt}$ (AB)", fontsize=10, rotation=270, labelpad=15.0)
    y_f = lambda y_i: y_i-vals.dm+2.5*np.log10(1+vals.z)
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.tick_params(axis='y', labelsize=10)
    ax2.set_yticks([-21,-20,-19,-18])
    ax2.set_yticklabels([-21,-20,-19,-18])
    ax2.plot([],[])

    # For the legend
    ax.plot([100,100],[150,150],ls='-',c='grey',label='SN2011kl')
    ax.plot([100,100],[150,150],ls='--',c=vals.gc,label='AT2018cow $g$',lw=0.5)
    ax.plot([100,100],[150,150],ls='--',c=vals.rc,label='AT2018cow $r$',lw=0.5)
    ax.plot([100,100],[150,150],ls='-',c=vals.gc,label='AT2020xnd $g$',lw=0.5)
    ax.plot([100,100],[150,150],ls='-.',c=vals.gc,label='AT2020mrf $g$',lw=1)
    #ax.plot([100,100],[150,150],ls='-',c=vals.rc,label='AT2020xnd $r$',lw=0.5)
    ax.plot([100,100],[150,150],ls=':',c=vals.rc,label='SN1998bw')
    ax.scatter(0,0,marker='s',c=vals.gc,edgecolor='k',label='AT2022tsd $g$')
    ax.scatter(0,0,marker='o',c=vals.rc,edgecolor='k',label='AT2022tsd $r$')
    ax.scatter(0,0,marker='D',c=vals.ic,edgecolor='k',label='AT2022tsd $i$')
    ax.scatter(0,0,marker='>',c=vals.wc,edgecolor='k',label='AT2022tsd $w$')
    ax.scatter(0,0,marker='<',c=vals.wc,edgecolor='k',label='AT2022tsd $o$')
    fig.legend(fontsize=8, bbox_to_anchor=(0.5,1.05), loc='upper center',
               ncol=6, handletextpad=0.1, columnspacing=0.8)

    # Formatting of the whole fig
    fig.text(0.5, 0, r'$\Delta t_\mathrm{obs}$ (days)', ha='center')
    axarr[0].spines['right'].set_visible(False)
    axarr[1].spines['left'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    axarr[0].tick_params(labelright='off')
    axarr[1].yaxis.tick_right()
    #axarr[1].set_yticks([])
    ax2.yaxis.tick_right()

    fig.subplots_adjust(wspace=0)
    plt.show()
    #plt.savefig("opt_lc.png", dpi=300, bbox_inches='tight', pad_inches=0.05)
    #plt.close()
