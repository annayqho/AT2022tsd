""" Optical rise time vs. peak luminosity
with a Nickel-56 curve """
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.cosmology import Planck18
from get_color_code import get_cc


ec,fc,msize,shape = get_cc()

iacol = '#003f5c'
cccol = '#58508d'
aftcol = '#bc5090'
llgrbcol = '#ff6361'
cowcol = '#ffa600'


def calc_Lpeak(Mni, tpeak):
    """
    Peak luminosity as a function of nickel mass and rise time

    Mni: nickel mass in units of Msun
    """
    tauni = 8.8
    tauco = 113.6
    Lpeak = 2E43 * (Mni) * (3.9*np.exp(-tpeak/tauni) + \
          0.678*(np.exp(-tpeak/tauco)-np.exp(-tpeak/tauni)))
    return Lpeak


def calc_Mej(tpeak):
    """
    The ejecta mass corresponding to rise time (Eq 2 of Kasen et al. 2017)
    We're assuming kappa = 0.1, V9=1

    Mej: in units of Msun
    tpeak: in units of days
    """
    Mej = (tpeak/14.5)**2
    return Mej


def plot_ztf(ax, background=False, shrink=1, text=True):
    """ Plot the ZTF sample """

    # Read table
    a = pd.read_csv("basic_info.csv")
    b = pd.read_csv("timescales.txt")

    # Get basic info
    names = a['Name'].values
    z = a['Redshift'].values
    cl = a['Class'].values
    m = a['mgpeak'].values
    em = a['emgpeak'].values

    # Get timescales
    name_key = b['name'].values
    grise = np.array([b['grise'].values[name_key==val][0] for val in names])

    # Munge the rises of the range events to have error bars
    tofix = np.array(['-' in val for val in grise])
    min_val = np.array([val.split('-')[0] for val in grise[tofix]]).astype(float)
    max_val = np.array([val.split('-')[1] for val in grise[tofix]]).astype(float)
    avgval = np.round((min_val+max_val)/2,2)
    errval = avgval-min_val
    repval = np.array(['%spm%s' %(i,j) for i,j in zip(avgval,errval)])
    grise[tofix] = repval

    # Munge the rises of the limit events
    tofix = np.array(['<' in val for val in grise])
    limval = np.array([val[1:] for val in grise[tofix]]).astype(float)
    grise[tofix] = np.array(['%spm%s' %(val,0) for val in limval])

    # Only keep events with known redshift
    keep = ~np.isnan(z)
    names = names[keep]
    z = z[keep].astype(float)
    cl = cl[keep]
    trise = np.array([float(val.split('pm')[0]) for val in grise[keep]])
    etrise = np.array([float(val.split('pm')[1]) for val in grise[keep]])

    M = m[keep].astype(float)-Planck18.distmod(z=z).value
    eM = em[keep].astype(float)

    i = 0
    for clname in np.unique(cl):
        if clname!='Unknown':
            choose = cl == clname

            # Scale by redshift
            x = trise[choose]/(1+z[choose])
            ex = etrise[choose]/(1+z[choose])

            y = M[choose]
            ey = eM[choose]

            zorder = 10
            if clname=='?':
                label = 'unclassified'
                zorder=0
            if np.logical_or.reduce((clname=='II', clname=='IIb', clname=='Ic', 
                clname=='Ib', clname=='Ibn', clname=='Ic-BL', clname=='IIn/Ibn')):
                clname = 'II'
                label = 'Core-collapse SN'
            else:
                label=clname

            # Go onto plotting
            if background is False:
                if clname=='AT2018cow-like':
                    col = cowcol
                else:
                    col = cccol
                for j,name in enumerate(names[choose]):
                    print(name)
                    ax.errorbar(
                            x[j], y[j], xerr=ex[j], yerr=ey[j],
                            label=None, mfc=col, mec=col,
                            c=col, fmt=shape[clname],
                            ms=msize[clname]/shrink, zorder=zorder)

    c = cowcol
    if background:
        c = 'lightgrey'
    # Plot label
    if text:
        ax.text(5, -22.5, "LFBOT", va='center', ha='right', color=c,
                fontweight='bold')
        xval = 10
        yval = -20.2
        xval = 12
        yval = -19.4


def plot_bts(ax):
    """ Plot the BTS sample with consistent coloring """
    dat = pd.read_csv("bts.csv")
    dur = dat['duration'].values
    Mpk = dat['peakabs'].values
    cl = dat['type'].values
    names = dat['ZTFID'].values

    keep = np.array([val!='-' for val in Mpk])
    dur = dur[keep]
    cl = cl[keep]
    names = names[keep]
    Mpk = Mpk[keep].astype(float)
    keep = np.array(['>' not in val for val in dur])
    dur = dur[keep].astype(float)
    Mpk = Mpk[keep]
    cl = cl[keep]
    names = names[keep]

    # Show them all in grey
    #ax.scatter(dur, Mpk, facecolor='grey', edgecolor='None', marker='s', alpha=0.3)

    ec,fc,msize,shape = get_cc()
    choose = np.logical_or.reduce((cl=='SN II', cl=='SN IIb', cl=='SN Ic', cl=='SN Ib',
            cl=='SN Ibn', cl=='SN Ic-BL'))
    ax.scatter(
            dur[choose], Mpk[choose],
            c=cccol, marker=shape['II'], zorder=2)
    ax.text(17,-15.7,'Core-collapse',c=cccol, fontweight='bold',ha='right')
    ax.text(13,-15.3,'SNe',c=cccol, fontweight='bold',ha='right')

    choose = cl=='SN Ia'
    ax.scatter(
            dur[choose], Mpk[choose],
            c=iacol, marker='.', zorder=0)
    ax.text(12,-20.2,'SN Ia',c=iacol, rotation=20, fontweight='bold')

    choose = ['SLSN' in val for val in cl]
    ax.scatter(
            dur[choose], Mpk[choose],
            c='grey', marker='x', zorder=0, s=10)
    ax.text(40,-22.5,'SLSN',c='grey', fontweight='bold')


if __name__=="__main__":
    # Initialize figure
    fig,ax = plt.subplots(1,1,figsize=(5,4))

    # Plot BTS sources
    plot_bts(ax)

    # Plot ZTF sources
    plot_ztf(ax, background=False, shrink=2, text=True)
    ax.set_ylabel("Peak absolute mag.", fontsize=14)
    ax.set_ylim(-15,-23)

    # Twin panel with peak lum
    ax2 = ax.twinx()
    # use the g-band effective wavelength: 4722.74 AA
    y_f = lambda y_i: 4E33*10**((y_i-4.77)/(-2.5))
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.plot([],[])
    ax2.set_yscale('log')
    ax2.tick_params(axis='both', labelsize=12)
    ax2.set_ylabel(
            r"Peak $\nu L_\nu$ (erg s$^{-1}$)", fontsize=14, rotation=270,
            labelpad=15.0)

    # Plot the Mni=Mej limit
    tpeak = np.linspace(1, 300)
    Mej = calc_Mej(tpeak)
    Lpeak = calc_Lpeak(Mej, tpeak)
    ax2.plot(tpeak, Lpeak, c='k', ls='--')

    # Formatting
    ax.set_xscale('log')
    ax.set_xlim(0.7,150)
    ax.set_xscale('log')
    ax.set_xlabel("Time above half-max (rest-frame days)", fontsize=14)
    ax.tick_params(axis='both', labelsize=12)

    # Save figure
    plt.tight_layout()
    plt.show()
