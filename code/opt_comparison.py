""" Compare the optical LC and spectrum to other classes of 
extragalactic transients.
Classes: 
    - CC SN
    - LFBOTs
    - TDEs
    - LGRBs 
    - LLGRBs 
"""


import matplotlib.pyplot as plt
from astropy.cosmology import Planck18
import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code)")
import vals
from get_opt import *
from helpers import bin_lc


def plot_22tsd(ax):
    """ Plot the extinction-corrected optical LC of AT2022tsd """
    # Get the LC of the baseline transient
    opt = get_full_opt()
    choose = np.logical_and(opt['isflare']==False, opt['mjdstart']<59856.4)
    opt_transient = opt[choose]
    dt = Time(opt_transient['mjdstart'].values, format='mjd').jd-vals.t0
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

        ax.errorbar(xb/(1+vals.z), yf-vals.dm, eyf, fmt=ms[i], color=cs[i])


def plot_18cow(ax):
    """ Plot the optical LC of AT2018cow """
    dat = pd.read_fwf("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/opt/at2018cow_photometry_table.dat")
    t0 = 58286
    x = Time(dat['MJD'].values.astype('float'), format='mjd').value
    dat['Emag'][dat['Emag']=='-']=99
    for i in np.arange(len(dat['ABMag'])):
        if '>' in dat['ABMag'][i]:
            dat['ABMag'][i] = dat['ABMag'][i][1:]
    M = dat['ABMag'].values.astype(float)-Planck18.distmod(z=0.0141).value
    cols = [vals.gc, vals.rc]
    for i,b in enumerate(['g', 'r']):
        choose = dat['Filt'].values==b
        ax.plot(x[choose]-t0, M[choose], c=cols[i], ls='--')
    ax.text(2, -18, '18cow', rotation=-35, fontsize=8)


def plot_20xnd(ax):
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
        y = mag[choose].astype(float)-Planck18.distmod(z=0.2442).value
        ax.plot(x/(1.2442), y, c=cols[i])
    ax.text(15, -17, '20xnd', rotation=-35, fontsize=8)


if __name__=="__main__":
    # Initialize a 2-panel figure
    fig,axarr = plt.subplots(1,2,figsize=(8,3))

    # Optical LC panel
    ax = axarr[0]
    plot_22tsd(ax)
    plot_18cow(ax)
    plot_20xnd(ax)

    ax.invert_yaxis()
    ax.set_xlim(-5, 40)
    ax.set_ylim(-15, -22)
    ax.set_xlabel("Rest-frame days")
    ax.set_ylabel("Absolute magnitude")

    plt.tight_layout()
    plt.show()
