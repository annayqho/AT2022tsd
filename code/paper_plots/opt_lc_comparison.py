""" Compare the optical LC and spectrum to other classes of transients.

Classes: 
    - CC SN
    - LFBOTs (LC done)
    - TDEs (22cmc LC is red, don't bother)
    - LGRBs (LCs are red, don't bother)
    - LLGRBs (98bw)
    - Stellar phenomena (CVs, novae)
"""


import matplotlib.pyplot as plt
from astropy.cosmology import Planck18
import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code)")
import vals
from get_opt import *
from helpers import *


def plot_22tsd(ax, show='absolute', offset=0):
    """ Plot the extinction-corrected optical LC of AT2022tsd """
    # Get the LC of the baseline transient
    opt = get_full_opt()
    choose = np.logical_and(opt['isflare']==False, opt['mjdstart']<59856.4)
    opt_transient = opt[choose]
    dt = Time(opt_transient['mjdstart'].values, format='mjd').jd-vals.t0+offset
    mag = opt_transient[choose]['mag_extcorr'].values
    emag = opt_transient[choose]['emag'].values
    filt = opt_transient['flt'].values

    # Add the extinction-corrected Keck data
    t = 59871.44599
    dtval = Time(t,format='mjd').jd-vals.t0

    g = 1.85
    eg = 0.14
    mg = -2.5*np.log10(g*1E-6)+8.90-vals.ext['g']
    emg = (2.5/np.log(10)) * (eg/g)
    dt = np.append(dt, dtval)
    mag = np.append(mag, mg)
    emag = np.append(emag, emg)
    filt = np.append(filt, 'g')

    i = 2.75
    ei = 0.15
    mi = -2.5*np.log10(i*1E-6)+8.90-vals.ext['i']
    emi = (2.5/np.log(10)) * (ei/i)
    dt = np.append(dt, dtval)
    mag = np.append(mag, mi)
    emag = np.append(emag, emi)
    filt = np.append(filt, 'i')

    # Keep the detections
    isdet = mag < 99
    dt = dt[isdet]
    mag = mag[isdet]
    emag = emag[isdet]
    filt = filt[isdet]

    # Plot the LC binned by one day
    ms = ['s', 'o', 'D']
    cs = [vals.gc, vals.rc, vals.ic]

    for i,b in enumerate(np.array(['g', 'r', 'i'])):
        # Get the LC in that filter
        choose = filt==b
        x = dt[choose]
        y = mag[choose]
        ey = emag[choose]

        # Convert to flux units
        flux = 3631*10**(y/(-2.5))
        eflux = (np.log(10)/2.5)*(ey)*flux

        # Bin 
        xb, yb, eyb = bin_lc(x, flux, eflux, 0.6)

        # Convert back to mag units
        yf = -2.5*np.log10(yb/3631)
        eyf = (eyb/yb) * (2.5/np.log(10))

        if show=='absolute':
            ax.errorbar(xb/(1+vals.z), yf-vals.dm, eyf, fmt=ms[i], color=cs[i],
                        zorder=10)
        elif show=='apparent':
            ax.errorbar(xb, yf, eyf, fmt=ms[i], color=cs[i], zorder=10)


def plot_18cow(ax, show='absolute', offset=0):
    """ Plot the optical LC of AT2018cow """
    dat = pd.read_fwf("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/opt/at2018cow_photometry_table.dat")
    t0 = 58286
    x = Time(dat['MJD'].values.astype('float'), format='mjd').value
    dat['Emag'][dat['Emag']=='-']=99
    for i in np.arange(len(dat['ABMag'])):
        if '>' in dat['ABMag'][i]:
            dat['ABMag'][i] = dat['ABMag'][i][1:]
    M = dat['ABMag'].values.astype(float)
    if show=='absolute':
        M = dat['ABMag'].values.astype(float)-Planck18.distmod(z=0.0141).value
    cols = [vals.gc, vals.rc]
    for i,b in enumerate(['g', 'r']):
        choose = dat['Filt'].values==b
        ax.plot(x[choose]-t0-1, M[choose]+offset, c=cols[i], ls='--', lw=0.7)
    #ax.text(2, -18, '18cow', rotation=-35, fontsize=8)


def plot_20xnd(ax, show='absolute', offset=0):
    """ Plot the extinction-corrected light curve of AT2020xnd """
    dat = pd.read_fwf("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/opt/at2020xnd_photometry_table.dat")
    mjd = dat['#MJD']
    mag = dat['mABc']
    emag = dat['merr']
    filt = dat['filt']
    t0 = 59135
    cols = [vals.gc, vals.rc]
    for i,b in enumerate(['g', 'r']):
        choose = np.logical_and(dat['filt'].values==b, dat['merr']!='---')
        x = mjd[choose]-t0
        y = mag[choose].astype(float)
        if show=='absolute':
            y = mag[choose].astype(float)-Planck18.distmod(z=0.2442).value+offset
        ax.plot(x[::2]/(1.2442), y[::2], c=cols[i], lw=0.5)
    #ax.text(15, -17, '20xnd', rotation=-35, fontsize=8)


def plot_98bw(ax, show='absolute', offset=0):
    """ Plot the LC of SN1998bw """
    dat = pd.read_csv(
            "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/opt/sn1998bw.dat",
            delimiter=';')
    dm = Planck18.distmod(z=0.0085).value
    jd = dat['JD'].values
    rband = dat['Rcmag'].values
    erband = dat['e_Rcmag'].values
    choose = rband!="     "
    jd = jd[choose]
    rband = rband[choose].astype(float)
    erband = erband[choose].astype(float)
    # Extinction is 0.127 in R-band in this direction
    y = rband
    if show=='absolute':
        y = rband-dm
    ax.plot(jd-jd[0], y-0.127+offset, color=vals.rc, lw=1, ls=':')
    #ax.text(30, -18.7, '98bw', rotation=-15, fontsize=8)


def plot_sn2011kl(ax, show='absolute', offset=0):
    """ Plot the LC of the ULGRB counterpart """
    z = 0.677
    dat = np.loadtxt("../2011kl.txt")
    t = dat[:,0]
    m = dat[:,1]
    dm = Planck18.distmod(z=z).value
    if show=='apparent':
        m = m+dm
    ax.plot(t/(1.677), m-offset, c='grey', lw=1, ls='-.')


def plot_cvs(ax):
    dat = pd.read_csv("bts_cv.txt")

    # Get rid of things with bad measurements
    rise = dat['rise'].values
    discard = np.array(['>' in val for val in rise])
    dat = dat[~discard]
    fade = dat['fade'].values
    discard = np.array(['>' in val for val in fade])
    dat = dat[~discard]

    names = dat['ZTFID'].values
    rise = dat['rise'].values.astype(float)
    fade = dat['fade'].values.astype(float)
    choose = np.logical_and(rise<10, fade<20)

    names = names[choose]

    s = logon()
    for i,name in enumerate(names):
        isalert,jd,mag,emag,filt,program,limjds,limmags,limfilts,limprogram =\
                get_lc(s, name)
        if max(jd)-min(jd)<50:
            choose = filt==1
            maxmag = min(mag)
            shift = 19.1-maxmag
            ax.plot(
                    jd[choose]-jd[choose][0], mag[choose]+shift, 
                    c='lightgrey', alpha=0.5, zorder=0)


def plot_novae(ax):
    dat = pd.read_csv("bts_novae.txt")

    # Get rid of things with bad measurements
    rise = dat['rise'].values
    discard = np.array(['>' in val for val in rise])
    dat = dat[~discard]
    fade = dat['fade'].values
    discard = np.array(['>' in val for val in fade])
    dat = dat[~discard]

    names = dat['ZTFID'].values
    rise = dat['rise'].values.astype(float)
    fade = dat['fade'].values.astype(float)
    choose = np.logical_and(rise<10, fade<20)

    names = names[choose]

    s = logon()
    for i,name in enumerate(names):
        isalert,jd,mag,emag,filt,program,limjds,limmags,limfilts,limprogram =\
                get_lc(s, name)
        if max(jd)-min(jd)<50:
            choose = filt==1
            maxmag = min(mag)
            shift = 19.1-maxmag
            ax.plot(
                    jd[choose]-jd[choose][0], mag[choose]+shift, 
                    c='lightgrey', alpha=0.5, zorder=0)


def compare_opt_lc():
    """ Three panels, comparing the optical LC to various transients """
    # Initialize a 3-panel figure
    fig,axarr = plt.subplots(1,3,figsize=(10,3))

    # LFBOTs
    ax = axarr[0]
    plot_22tsd(ax)
    plot_18cow(ax)
    plot_20xnd(ax)
    plot_98bw(ax)
    plot_sn2011kl(ax)
    ax.text(0.95, 0.95, "Extragalactic", va='top',
            ha='right', transform=ax.transAxes, fontsize=9)

    # Format the panel
    ax.invert_yaxis()
    ax.set_xlim(-5, 40)
    ax.set_ylim(-15, -22)
    ax.set_xlabel("Rest-frame days since explosion")
    ax.set_ylabel("Absolute magnitude")

    # Galactic phenomena
    ax = axarr[1]
    plot_22tsd(ax, show='apparent')
    ax.invert_yaxis()
    ax.set_xlabel("Days")
    ax.set_ylabel("Apparent magnitude")
    plot_cvs(ax)
    ax.text(0.95, 0.95, "CVs / Dwarf Novae", va='top',
            ha='right', transform=ax.transAxes, fontsize=9)
    ax.set_ylim(22.5,18.9)

    ax = axarr[2]
    plot_22tsd(ax, show='apparent')
    ax.invert_yaxis()
    ax.set_xlabel("Days")
    plot_novae(ax)
    ax.text(0.95, 0.95, "Novae", va='top',
            ha='right', transform=ax.transAxes, fontsize=9)

    plt.tight_layout()
    #plt.show()
    plt.savefig("compare_opt_lc.png", dpi=200)
    plt.close()

