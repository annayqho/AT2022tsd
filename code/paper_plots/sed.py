""" Plot the X-ray to radio SED
Currently the only good date for this is 27-28d.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from astropy.time import Time
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
import vals
from get_opt import *
from get_xray import load_swift
from get_radio_at2022tsd import get_radio


def plot_xray(ax, flux):
    """ Take the integrated flux across 0.3-10 keV,
    use that and the geometric mean of the frequency and the spectral index
    to solve for the normalization coefficient for your spectrum

    Parameters
    ----------
    y: flux

    2.4E18 Hz corresponds to 10 keV
    7.2E16 Hz corresponds to 0.3 keV
    """
    nu0 = np.sqrt(0.3*10) * (2.4E18/10)
    alpha = 1.01
    nu2 = 2.4E18
    nu1 = 7.2E16
    A = flux*(1-alpha) / (nu0**alpha * (nu2**(1-alpha) - nu1**(1-alpha)))
    print(A)
    xplot = np.linspace(7.2E16, 2.4E18, 1000)
    yplot = A*(xplot/nu0)**(-alpha)
    #print(xplot, xplot*yplot)
    ax.plot(xplot, xplot*yplot, c='k', ls='-')


def quiescent_sed():
    # Initialize
    fig,ax = plt.subplots(1,1,figsize=(4,3))

    # Add the radio to the plot
    dat = get_radio()
    dat['dt'] = Time(dat['Date'].values, format='isot').jd-vals.t0
    choose = np.logical_and(dt>27, dt<28) # 27.14 VLA, 27.70 NOEMA
    y = dat['Flux'][choose]*1E-3*1E-23*4*np.pi*(vals.dL_cm)**2
    ey = dat['eFlux'][choose]*1E-3*1E-23*4*np.pi*(vals.dL_cm)**2
    ax.errorbar(
            dat['Freq_Obs'][choose]*1E9, y*dat['Freq_Obs'][choose]*1E9, 
            ey*dat['Freq_Obs'][choose]*1E9, fmt='o', c='k')

    # Add the optical to the plot 
    dat = get_full_opt()
    dat['dt'] = (Time(dat['mjdstart'], format='mjd').jd-vals.t0)
    dt = dat['dt']
    choose_dt = np.logical_and(dt>27, dt<28)
    # use the NOT gri data
    choose_det = np.logical_and(dat['emag']<99, dat['#instrument']=='NOT/ALFOSC')
    choose = np.logical_and(choose_dt, choose_det)
    freq = np.zeros(len(dat[choose]['flt']))
    freq[dat[choose]['flt']=='i'] = 3E18/vals.sdss_pivot['i']
    freq[dat[choose]['flt']=='r'] = 3E18/vals.sdss_pivot['r']
    freq[dat[choose]['flt']=='g'] = 3E18/vals.sdss_pivot['g']
    y = dat['flux'][choose]*1E-6*1E-23*4*np.pi*(vals.dL_cm)**2
    ey = dat['unc'][choose]*1E-6*1E-23*4*np.pi*(vals.dL_cm)**2
    ax.errorbar(freq, y*freq, ey*freq, fmt='o', c='k', lw=0.2)

    # Add the X-ray to the plot
    df = load_swift()
    t = Time(df['!MJD    '].values, format='mjd')
    dt = t.jd-vals.t0
    et = (df['T_+ve   '].values)
    L = df['L'].values
    eL = df['Lpos'].values
    choose = np.logical_and(dt>27, dt<28)
    plot_xray(ax, L[choose][0])

    # Add a line
    xvals = np.linspace(1E10,1E13)
    yvals = (4E40)*(xvals/1E11)**2
    plt.plot(xvals, yvals, c='lightblue', label=r'$\nu L_\nu \propto \nu^2$')

    # Add a line
    xvals = np.linspace(1E10,1E17)
    yvals = (3E42)*(xvals/4.6E14)**0.3
    plt.plot(xvals, yvals, c='lightblue', 
             ls='--', label=r'$\nu L_\nu \propto \nu^{0.3}$')

    ax.text(0.05, 0.95, r'$\Delta t_\mathrm{rest}\approx30$d', 
            transform=ax.transAxes, ha='left', va='top')

    # Formatting
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(7E9, 4E18)
    ax.set_ylim(5E38,2E43)
    ax.legend(loc='lower right')

    ax.set_xlabel(
            r"$\nu_\mathrm{obs}$",fontsize=11,
            fontname='sans-serif')

    ax.set_ylabel(r"$\nu L_\nu$", fontsize=11,
            fontname='sans-serif')

    ax.tick_params(axis='both', labelsize=11)
    #plt.tight_layout()
    #plt.show()

    plt.savefig("sed.png", dpi=200, bbox_inches='tight', pad_inches=0.1)
    plt.close()
