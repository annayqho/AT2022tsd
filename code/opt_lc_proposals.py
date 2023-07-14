""" Plot the optical light curve """
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.cosmology import Planck15
from read_ipac_forced_phot import get_lc


def plot_lc(ax,f,t,m,em,markers=False,lines=False,lw=1):
    cols = ['Crimson', 'Aquamarine', 'Goldenrod']
    symbols = ['o', 's', '>']
    for i,filt in enumerate(np.array(['r', 'g', 'i'])):
        choose = np.logical_and(em<99, f==filt)
        if markers:
            ax.errorbar(t[choose], m[choose], em[choose], 
                         fmt=symbols[i], c=cols[i])
        if lines:
            ax.plot(t[choose], m[choose], c=cols[i], lw=lw)
        first_det_t = t[choose][0]
        first_det_m = m[choose][0]
        choose = np.logical_and(em==99, f==filt)
        if sum(choose)>0:
            if markers:
                ax.scatter(t[choose], m[choose], marker='v', c=cols[i])
            last_lim_t = t[choose][-1]
            last_lim_m = m[choose][-1]
            if lines:
                ax.plot([last_lim_t, first_det_t], [last_lim_m, first_det_m],
                         ls='--', c=cols[i], lw=lw)


def lc_AT2018cow(ax):
    dat = pd.read_fwf("../data/at2018cow_photometry_table.dat")
    t0 = 58287.1500
    x = Time(dat['MJD'].values.astype('float'), format='mjd').value
    dat['Emag'][dat['Emag']=='-']=99
    for i in np.arange(len(dat['ABMag'])):
        if '>' in dat['ABMag'][i]:
            dat['ABMag'][i] = dat['ABMag'][i][1:]
    M = dat['ABMag'].values.astype(float)-Planck15.distmod(z=0.0141).value
    plot_lc(ax,
            dat['Filt'].values, x-t0, M, dat['Emag'].values.astype(float), 
            lines=True, lw=0.5)
    ax.text(6, -17, 'AT2018cow')


def lc_AT2022tsd(ax):
    dm = Planck15.distmod(z=0.256).value

    # LC of our object
    t0 = Time("2022-09-07", format='isot')
    dat = pd.read_csv("../data/opt/opt.txt")

    # Plot the baseline LC
    x = Time(dat['Date'].values.astype(float), format='mjd')-t0
    choose = x < 40
    plot_lc(ax,
            dat['Filt'].values[choose], x.value[choose], 
            dat['Mag'].values[choose], dat['Emag'].values[choose], 
            markers=True, lines=True, lw=2)

    # Plot the flares
    choose = x > 40
    x = Time(dat['Date'].values.astype(float), format='mjd')-t0
    ax.errorbar(x.value[choose], 
               dat['Mag'].values[choose],
               dat['Emag'].values[choose], fmt='s', c='Aquamarine')
    #ax.text(12,-19.2,'AT2022tsd')

    # Get the IPAC LC
    jd,filt,mag,emag = get_lc()
    gdet = np.logical_and(filt=='g', emag<99)
    rdet = np.logical_and(filt=='r', emag<99)
    idet = np.logical_and(filt=='i', emag<99)
    glim = np.logical_and(filt=='g', emag==99)
    rlim = np.logical_and(filt=='r', emag==99)
    ilim = np.logical_and(filt=='i', emag==99)
    ax.errorbar(
            jd[gdet]-t0.jd, mag[gdet], emag[gdet], fmt='s', c='Aquamarine') 
    ax.errorbar(
            jd[rdet]-t0.jd, mag[rdet], emag[rdet], fmt='o', c='Crimson')
    ax.errorbar(
            jd[idet]-t0.jd, mag[idet], emag[idet], fmt='>', c='Goldenrod')
    ax.scatter(jd[glim]-t0.jd, mag[glim], marker='v', facecolor='white',
               edgecolor='Aquamarine',lw=0.5)
    ax.scatter(jd[rlim]-t0.jd, mag[rlim], marker='v', edgecolor='Crimson',
               facecolor='white',lw=0.2)
    ax.scatter(jd[ilim]-t0.jd, mag[ilim], marker='v', edgecolor='Goldenrod',
               facecolor='white',lw=0.3)


def lc_flares():
    fig,axarr = plt.subplots(1,2,figsize=(6,2.5), sharey=True)
    dm = Planck15.distmod(z=0.256).value

    # LC of our object
    t0 = Time("2022-09-07", format='isot')
    dat = pd.read_csv("../data/opt.txt")
    tel = dat['Tel'].values
    x = Time(dat['Date'].values.astype(float), format='mjd')-t0

    # Plot the flares
    ax = axarr[0]
    choose = np.logical_and(x > 40, tel=='IMACS')
    ax.errorbar((x.value[choose]-x.value[choose][0])*24*60, 
               dat['Mag'].values[choose],
               dat['Emag'].values[choose], fmt='s', c='Aquamarine')
    ax.set_ylabel("Apparent Mag.")

    ax = axarr[1]
    choose = np.logical_and(x > 40, tel=='LT')
    ax.errorbar((x.value[choose]-x.value[choose][0])*24*60, 
               dat['Mag'].values[choose],
               dat['Emag'].values[choose], fmt='s', c='Aquamarine')
    ax.invert_yaxis()
    # Second axis
    ax2 = ax.twinx()
    y_f = lambda y_i: y_i-dm
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.set_ylabel("Absolute Mag.", rotation=270, labelpad=15.0)
    ax2.plot([],[])

    for ax in axarr:
        ax.set_xlabel("Minutes")

    plt.tight_layout()
    #plt.show()
    plt.savefig("lc_flares.png", dpi=200)


def lc_flares_nulnu():
    fig,ax = plt.subplots(1,1,figsize=(3.5,2.5), sharey=True)
    dm = Planck15.distmod(z=0.256).value

    # LC of our object
    t0 = Time("2022-09-07", format='isot')
    dat = pd.read_csv("../data/opt.txt")
    tel = dat['Tel'].values
    x = Time(dat['Date'].values.astype(float), format='mjd')-t0

    # Plot the flares
    choose = np.logical_and(x > 40, tel=='IMACS')
    mag = dat['Mag'].values[choose]
    Mag = mag-dm
    lum = 4E33*10**((Mag-4.77)/(-2.5))

    ax.scatter((x.value[choose]-x.value[choose][0])*24*60, 
                lum, marker='s', c='Aquamarine')
    ax.set_ylabel(r"$\nu L_\nu$ (erg/s)")
    ax.set_yscale('log')

    ax.set_xlabel("Minutes")

    plt.tight_layout()
    #plt.show()
    plt.savefig("lc_flares_lum.png", dpi=200)
    plt.close()



def lc_ibc(ax):
    """ LC of a Type Ibc SN """
    dat = np.loadtxt("../data/ibc_template_R.txt", delimiter=',')
    y = dat[:,1]
    ax.plot(dat[:,0], y-19, color='grey', lw=1, ls=':', label='Ibc Template')



def plot_main_lc():
    dm = Planck15.distmod(z=0.256).value
    #fig,ax = plt.subplots(1,1,figsize=(6,2.5))
    fig,ax = plt.subplots(1,1,figsize=(5,2.5))
    lc_AT2022tsd(ax)
    #lc_AT2018cow(ax)
    #lc_ibc(ax)
    #plt.legend(loc='upper right')
    plt.xlabel("Days Since Peak")
    plt.ylabel("Absolute Mag")
    #ax.set_xlim(-10,105)
    ax.set_xlim(26.6,29.6)
    ax.set_ylim(25,19)

    # Second axis
    ax2 = ax.twinx()
    y_f = lambda y_i: y_i-dm
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.set_ylabel("Absolute Mag.", rotation=270, labelpad=15.0)
    ax2.plot([],[])

    plt.tight_layout()
    plt.show()
    #plt.savefig("opt_lc.png", dpi=200)
    #plt.close()


if __name__=="__main__":
    #lc_flares_nulnu()
    plot_main_lc()
