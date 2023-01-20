""" 
Plot of luminosity over time
"""


import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
from astropy.cosmology import Planck15
from astropy.time import Time
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
from get_radio import *
from scale_fluxes import sma_lc
from read_table import *
from format import *

form = get_format()

    

def sn1993J(ax, col, legend):
    """ SN 1993J from Weiler et al. 
    This is the peak of the 99.4 GHz light curve
    There is also a 110 GHz point, but only one,
    so I can't get the peak.
    1.4 GHz / 20 cm, 4.9 GHz / 6 cm
    values come from the best-fit model,
    but by eye they are clearly pretty close
    """
    freq = 99.4E9
    d = 1.1E25
    nu, dt, f, ef, islim = read_1993J_high_freq()
    l = f*1E-3*1E-23*4*np.pi*d**2
    choose = np.logical_and(~islim, nu==freq)
    ax.plot(dt[choose], l[choose], c=col, label=None, ls='--')
    ax.scatter(dt[choose], l[choose], c=col, label=legend, marker='s')
    #ax.text(dt[choose][5], l[choose][5]*1.4, 'SN1993J', 
    #        color=col, va='bottom', ha='center', fontsize=form['font_small'])


def sn2011dh(ax, col, legend):
    """ SN 2011dh
    Horesh et al. 2013
    Krauss et al. 2012
    M51: d = 8.03 Mpc; expl date May 31.58
    """
    d = 2.5E25

    # HIGH FREQUENCY
    # use two freq: 107E9 and 93E9
    dt, nu, f, ef, islim = read_2011dh()
    choose = np.logical_and(~islim, np.logical_or(nu==107E9, nu==93E9))
    l = f*1E-3*1E-23*4*np.pi*d**2
    ax.plot(dt[choose], l[choose], c=col, ls='--')
    ax.scatter(dt[choose], l[choose], c=col, marker='s')
    #ax.text(dt[choose][0]/1.1, l[choose][0], 'SN2011dh', c=col, ha='right',
    #        fontsize=form['font_small'])


def ptf11qcj(ax, col, legend):
    """ PTF 11qcj
    Corsi et al. 2014
    """
    d = 124 * 3.086E24

    # HIGH FREQUENCY (Carma)
    t0 = 55857.543-15
    dt = [55884.803-t0, 55891.640-t0]
    f = np.array([3.96, 3.62])
    ef = [0.88, 0.75]
    nu = 93E9

    l = f*1E-3*1E-23*4*np.pi*d**2
    ax.plot(dt, l/1.2, c=col, ls='--', zorder=10)
    ax.scatter(dt, l/1.2, c=col, marker='s', zorder=10)
    #ax.text(dt[0]/1.2, l[0]/1.5, 'PTF11qcj', color=col, va='top', ha='center',
            #fontsize=form['font_small'])


def sn2008d(ax, col, legend):
    """ SN 2008D
    Soderberg et al. 
    """
    d = 29.9 * 3.086E24

    # HIGH FREQUENCY (Carma)
    dt = [4.94, 6.84]
    f = np.array([3.2, 0.6])
    ef = [0.7, 0.3]
    nu = 95E9

    l = f*1E-3*1E-23*4*np.pi*d**2
    ax.scatter(dt, l/1.2, c=col, marker='s')
    ax.plot(dt, l/1.2, c=col, ls='--')
    #ax.text(dt[0]/1.1, l[0]/1.2, 'SN2008D', color=col, ha='right',
    #        fontsize=form['font_small'])


def sn2020oi(ax, col, legend):
    """ Maeda et al. 2021 """
    d = 15.5 * 3.086E24
    dt = np.array([5.4, 8.4, 18.3, 51.3])
    fnu = np.array([1.3, 1.22, 0.196, 0.115])
    l = fnu*1E-3*1E-23*4*np.pi*d**2
    ax.scatter(dt, l, color=col, marker='s')
    #ax.text(dt[-1]*1.1, l[-1], 'SN2020oi', color=col, ha='left',
    #        fontsize=form['font_small'])
    ax.plot(dt, l, c=col, ls='--')


def grb030329(ax, col, legend):
    """ Sheth et al. 2003 
    Berger 2003
    Van der Horst et al. 2008
    
    Explosion day was obviously 03/29
    """
    z = 0.1686
    d = Planck15.luminosity_distance(z=z).cgs.value

    # HIGH FREQUENCY

    # Sheth
    dat = np.loadtxt("030329_dat.txt", delimiter=',', dtype=str)
    freq = 100E9
    t = dat[:,0].astype(float)
    lum = dat[:,1].astype(float) * 1E-3 * 1E-23 * 4 * np.pi * d**2
    ax.plot(t, lum, c=col)
    #ax.text(t[1], lum[1]/1.1, 'GRB030329', 
    #        color=col, va='top', ha='center', fontsize=form['font_small'])


def grb181201A(ax, col, legend):
    """ Laskar et al. 2018
    """
    z = 0.450
    d = Planck15.luminosity_distance(z=z).cgs.value

    freq = 97
    t = np.array([8.847e-01, 1.917e+00, 3.877e+00, 8.661e+00, 2.981e+01])
    flux = np.array([3.413, 1.987, 1.199, 0.624, 0.259])
    lum = flux * 1E-3 * 1E-23 * 4* np.pi*d**2
    ax.plot(t, lum, c=col)
    #ax.text(t[-1], lum[-1]*1.1, 'GRB181201A', 
    #        color=col, va='bottom', ha='left', fontsize=form['font_small'])


def grb161219B(ax, col, legend):
    """ Laskar et al. """
    z = 0.1475
    d = Planck15.luminosity_distance(z=z).cgs.value

    freq = 100
    t = np.array([1.30, 3.30, 8.31, 24.45, 78.18])
    flux = np.array([1244, 897, 500, 285, 51])
    lum = flux * 1E-6 * 1E-23 * 4* np.pi*d**2
    ax.plot(t, lum, c=col)
    #ax.text(t[-1], lum[-1]/1.1, 'GRB161219B', 
    #        color=col, va='top', ha='center', fontsize=form['font_small'])


def grb130427A(ax, col, legend):
    """ Perley et al
    They have data from CARMA/PdBI at 90 GHz (3mm)
    and also CARMA at 85.00 GHz, which is close enough
    But by the time they caught it, it was fading
    """
    z = 0.340
    d = Planck15.luminosity_distance(z=z).cgs.value

    freq = 93E9
    obs_t_1 = np.array([0.77, 1, 1.91, 2.8]) 
    obs_flux_1 = np.array([3416, 2470, 1189, 807]) * 1E-6
    obs_lum_1 = obs_flux_1 * 1E-23 * 4 * np.pi * d**2

    freq = 85E9 
    obs_t_2 = np.array([0.81, 3.58, 6.41, 10.36, 23.52])
    obs_flux_2 = np.array([3000, 903, 588, 368, 197]) * 1E-6
    obs_lum_2 = obs_flux_2 * 1E-23 * 4 * np.pi * d**2

    t = np.hstack((obs_t_1, obs_t_2))
    lum = np.hstack((obs_lum_1, obs_lum_2))
    order =np.argsort(t)

    ax.plot(obs_t_2, obs_lum_2, c=col, label=legend)
    #ax.text(t[2], lum[2]/1.6, 'GRB130427A', 
    #        color=col, va='top', ha='center', fontsize=form['font_small'])


def j1644(ax, col, legend):
    """ Zauderer
    They have data from CARMA at 94 GHz (3mm)
    and also CARMA at 87.00 GHz, which is close enough
    But by the time they caught it, it was fading
    """
    z = 0.354
    d = Planck15.luminosity_distance(z=z).cgs.value

    freq = 94E9
    obs_t_1 = np.array([1.85, 19.06])
    obs_flux_1 = np.array([15.7, 10.7]) * 1E-3
    obs_lum_1 = obs_flux_1 * 1E-23 * 4 * np.pi * d**2

    freq = 87E9 
    obs_t_2 = np.array([5.14, 6.09, 7.18, 9.09, 14.61, 19.06, 22.07])
    obs_flux_2 = np.array([18.6, 21.7, 14.6, 15.1, 10.4, 9.36, 5.49])*1E-3
    obs_lum_2 = obs_flux_2 * 1E-23 * 4 * np.pi * d**2

    t = np.hstack((obs_t_1, obs_t_2))
    lum = np.hstack((obs_lum_1, obs_lum_2))
    order =np.argsort(t)

    ax.scatter(t[order], lum[order], marker='o',
            facecolor='white', edgecolor=col, label=legend)
    ax.plot(t[order], lum[order], c=col, label=None, lw=1)
    #ax.text(t[order][-1]*1.2, lum[order][-1], 'J1644+57', 
    #        ha='left', va='bottom',
    #        fontsize=form['font_small'], color=col)


def igr(ax, col, legend):
    """ IGR J12580+0134
    They have data from Planck at 100 GHz
    Discovered by INTEGRAL (https://www.astronomerstelegram.org/?read=3108)

    First detection: 2011 Jan 2-11
    Last non-detection: 2010 Dec 30 to 2011 Jan 2

    So... the dt is something like 1 day to 12 days?
    """
    d = 17*3.086E24 # Mpc to cm

    t = 10 # estimate
    freq = 100E9
    lum = 640*1E-3*1E-23*4*np.pi*d**2

    ax.scatter(t, lum, marker='o',
            facecolor='white', edgecolor=col, label=legend)
    #ax.text(t/1.2, lum, 'IGR J12580', ha='right', va='center',
    #        fontsize=form['font_small'], color=col)


def sn1998bw(ax, col, legend):
    """ SN 1998bw
    
    This is a bit complicated because there are two peaks in the light curve,
    but I am using a frequency that has a main peak rather than a frequecy
    with two clear distinct peaks...
    """
    d = 1.17E26 # cm
    nu = 150E9
    t = np.array([12.4])
    f = np.array([39])
    lum = f*1E-3*1E-23*4*np.pi*d**2
    ax.scatter(t, lum, c=col, label=legend)
    #ax.text(t, lum*1.2, 'SN1998bw', color=col, va='bottom', ha='center',
    #        fontsize=form['font_small'])


def sn2017iuk(ax, col, legend):
    """ SN 2017iuk
    """
    d = Planck15.luminosity_distance(z=0.0368).cgs.value
    nu = 92E9 # Band 3
    t = np.array([6.10])
    f = np.array([28])
    l = f*1E-3*1E-23*4*np.pi*d**2
    ax.scatter(t,l,c=col)
    #ax.text(t/1.1, l, 'SN2017iuk', color=llgrb_col, ha='right',
    #        fontsize = form['font_small'])


def sn2020bvc(ax, col, legend):
    """ SN 2020bvc
    """
    d = Planck15.luminosity_distance(z=0.0252).cgs.value
    nu = 230E9 # Band 3
    t = 5.8
    f = 0.25
    l = f*1E-3*1E-23*4*np.pi*d**2
    ax.scatter(t,3*l,c=col)
    ax.arrow(t,3*l,0,-(3*l)/2,
            length_includes_head=True,head_length=(3*l)/6,
            head_width=t/7,color=col)
    #ax.text(t*1.1, l, 'SN2020bvc (230 GHz)', color=llgrb_col, ha='left')


def at2018cow(ax, col, legend):
    """ 230 GHz values """
    dcm = Planck15.luminosity_distance(z=0.0141).cgs.value
    a, b, c = sma_lc()
    dt, f, ef = b
    ef_comb = np.sqrt(ef**2 + (0.15*f)**2)
    nu = 231.5E9
    lum = f * 1E-3 * 1E-23 * 4 * np.pi * dcm**2
    ax.scatter(dt, lum, c=col, marker='D', label=legend, zorder=10)
    ax.plot(dt, lum, c=col, ls=':', lw=2, label=None, zorder=10)
    ax.text(
            dt[-1]*1.1, lum[-1]/1.5, 'AT2018cow', ha='center', va='top',
            color=col, fontsize=form['font_small'])



def at2020xnd(ax, col, legend):
    dcm = Planck15.luminosity_distance(z=0.2442).cgs.value
    dt = [17,24,31,39,46,67]
    fnu = np.array([305,648,1030,868,558,330])
    efnu = np.array([57,44,44,46,38,48])
    nu = np.array([94,94,94,94,94,79])
    lum = fnu * 1E-6 * 1E-23 * 4 * np.pi * dcm**2 
    ax.scatter(dt, lum, c=col, marker='D', label=legend, zorder=10)
    ax.plot(dt, lum, c=col, ls=':', lw=2, label=None, zorder=10)
    ax.text(
            dt[2]*1.4, lum[2]*1.3, 'AT2020xnd', ha='center', va='bottom',
            color=col, fontsize=form['font_small'])


def at2022tsd(ax, col, legend):
    dcm = Planck15.luminosity_distance(z=0.256).cgs.value
    dt = [27, 35, 51]
    fnu = np.array([0.245, 0.316, 0.299])
    efnu = np.array([0.065, 0.078, 0.093])
    nu = np.array([117, 117, 117])
    lum = fnu * 1E-3 * 1E-23 * 4 * np.pi * dcm**2 
    ax.scatter(
            dt, lum, facecolor=col, marker='D', 
            label=legend, zorder=10, edgecolor='k')
    ax.plot(dt, lum, c=col, ls=':', lw=2, label=None, zorder=10)
    ax.text(
            dt[-1]*2, lum[-1]/1.5, 'AT2022tsd', ha='center', va='top',
            color=col, fontsize=form['font_small'], fontweight='bold')


def at2022cmc(ax, col, legend):
    dcm = Planck15.luminosity_distance(z=1.193).cgs.value
    t0 = Time(59621.4, format='mjd')
    dt = Time(
            np.array(['2022-02-18','2022-02-20', '2022-02-22', '2022-02-24',
             '2022-02-25', '2022-02-26', '2022-02-28', '2022-03-02',
             '2022-03-06']), format='isot')-t0
    dt = dt.value

    fnu = np.array([9.1, 7.4, 7.9, 7.4, 6.6, 4.65, 5.0, 4.0, 3.0])
    efnu = np.array([0.9]*len(fnu))
    lnu = fnu*1E-3*1E-23*4*np.pi*dcm**2

    ax.scatter(dt, lnu, marker='o', edgecolor=col, facecolor='white')
    ax.plot(dt,lnu,c=col)
    #ax.text(
    #        dt[-1]*1.01, lnu[0], '22cmc (230 GHz)', fontsize=10, 
    #        ha='left', va='center', color='red', fontweight='bold')



def run(ax):
    props = dict(boxstyle='round', facecolor='white')

    sn_col = form['colors']['5'][0]
    llgrb_col = form['colors']['5'][1]
    cow_col = form['colors']['5'][2]
    lgrb_col = form['colors']['5'][3]
    tde_col = form['colors']['5'][4]

    # Category: TDEs
    j1644(ax, tde_col, legend='TDE')
    igr(ax, tde_col, legend=None)
    at2022cmc(ax, tde_col, legend=None)

    # First category: long-duration GRBs
    grb130427A(ax, lgrb_col, legend='LGRB')
    grb030329(ax, lgrb_col, None)
    grb181201A(ax, lgrb_col, None)
    grb161219B(ax, lgrb_col, None)

    # Second category: low-luminosity GRBs
    sn1998bw(ax, llgrb_col, legend='LLGRB')
    sn2017iuk(ax, llgrb_col, None)
    #sn2006aj(ax, llgrb_col, None)

    # Third category: Cow-like
    at2018cow(ax, cow_col, 'LFBOT')
    at2020xnd(ax, cow_col, None)
    at2022tsd(ax, cow_col, None)

    # Final category: 'ordinary' CC SNe
    sn1993J(ax, sn_col, 'CC SN')
    sn2011dh(ax, sn_col, None)
    sn2020oi(ax, sn_col, None)
    ptf11qcj(ax, sn_col, None)
    sn2008d(ax, sn_col, None)

    ax.legend(fontsize=8, loc='lower left')


if __name__=="__main__":
    fig,ax = plt.subplots(1,1,figsize=(3,4))
    run(ax)
    # Formatting
    ax.set_ylabel(
            r"$L_{\nu}$ (erg$\,$s$^{-1}$Hz$^{-1}$)", 
            fontsize=form['font_med'])
    ax.set_title(
            r"$\nu_\mathrm{obs} \gtrsim 100\,\mathrm{GHz}$", 
            fontsize=form['font_med'])
    ax.tick_params(axis='both', labelsize=form['font_small'])
    ax.set_xlim(0.7, 300) 
    ax.set_ylim(1E25, 2E33)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r"$\Delta t_\mathrm{obs}$ (d)", fontsize=form['font_med'])

    plt.tight_layout()
    plt.show()
    #plt.savefig(
    #        "mm_lc_100ghz.png", dpi=300, 
    #        bbox_inches='tight', pad_inches=0.1)
    #plt.close()
