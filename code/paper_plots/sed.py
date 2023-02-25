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
    y = xplot*yplot
    ax.plot(xplot, y, c='k', ls='-')
    ax.text(xplot[150], y[0]/1.2, 'Swift/XRT', fontsize=8, ha='center', va='top')


def quiescent_sed(ax):
    start_time = 26
    end_time = 28.7

    # The NOT optical data at dt=27.86512 (2022-10-04T02:10:46.848)
    # still appears blue. 
    # The LT optical data at 2022-10-04T02:44:11.328 (dt=27.88832) also 
    # appears blueish.
    # I think the flaring hasn't started yet.

    # The first flare is with ZTF, dt=28.186545 (2022-10-04T09:53:37.968)

    # Add the radio to the plot
    dat = get_radio()
    dt =  Time(dat['Date'].values.astype(str), format='isot').jd-vals.t0
    # There's only one radio detection before 28d...
    choose = np.logical_and(dt>start_time, dt<end_time) 
    x = dat['Freq_Obs'][choose].values *1E9
    y = dat['Flux'][choose].values*1E-3*1E-23*4*np.pi*(vals.dL_cm)**2 * x 
    ey = dat['eFlux'][choose].values*1E-3*1E-23*4*np.pi*(vals.dL_cm)**2 * x
    ax.errorbar(x, y, ey, fmt='o', c='k')
    ax.text(x[0]*1.5, y[0], 'VLA', ha='left', va='center', fontsize=8)
    ax.text(x[-1]*1.8, y[-1], 'NOEMA', ha='left', va='center', fontsize=8)

    # Add the optical to the plot 
    dat = get_full_opt()
    dat['dt'] = (Time(dat['mjdstart'], format='mjd').jd-vals.t0)
    dt = dat['dt']
    choose_dt = np.logical_and(dt>start_time, dt<end_time)
    # use the NOT gri data
    choose_det = np.logical_and(dat['emag']<99, dat['#instrument']=='NOT/ALFOSC')
    choose = np.logical_and(choose_dt, choose_det)
    freq = np.zeros(len(dat[choose]['flt']))
    freq[dat[choose]['flt']=='i'] = 3E18/vals.sdss_pivot['i']
    freq[dat[choose]['flt']=='r'] = 3E18/vals.sdss_pivot['r']
    freq[dat[choose]['flt']=='g'] = 3E18/vals.sdss_pivot['g']
    y = dat['flux'][choose].values*1E-6*1E-23*4*np.pi*(vals.dL_cm)**2 * freq
    # correct for Milky Way extinction
    fac = np.zeros(len(dat[choose]['flt']))
    fac[dat[choose]['flt']=='i'] = 10**(vals.ext['i']/2.5)
    fac[dat[choose]['flt']=='r'] = 10**(vals.ext['r']/2.5)
    fac[dat[choose]['flt']=='g'] = 10**(vals.ext['g']/2.5)
    y = y*fac
    ey = dat['unc'][choose].values*1E-6*1E-23*4*np.pi*(vals.dL_cm)**2 * freq * fac
    ax.errorbar(freq, y, ey, fmt='o', c='k', lw=0.2)
    ax.text(freq[-1]*1.8, y[-1], 'NOT', fontsize=8, ha='left', va='center')

    # Add the X-ray to the plot
    # The first Swift data is at dt=28 days
    df = load_swift()
    t = Time(df['!MJD    '].values, format='mjd')
    dt = t.jd-vals.t0
    et = (df['T_+ve   '].values)
    L = df['L'].values
    eL = df['Lpos'].values
    choose = np.logical_and(dt>start_time, dt<end_time)
    plot_xray(ax, L[choose][0])

    # Add a line
    xvals = np.linspace(1E10,1E13)
    yvals = (4.6E40)*(xvals/8.38E10)**2.5
    ax.plot(xvals, yvals, c='lightblue', label=r'$\nu L_\nu \propto \nu^{2.5}$')

    # Add a line
    xvals = np.linspace(1E10,1E17)
    yvals = (3E42)*(xvals/4.6E14)**0.5
    ax.plot(xvals, yvals, c='lightblue', 
             ls='--', label=r'$\nu L_\nu \propto \nu^{0.5}$')

    ax.text(0.05, 0.95, r'$\Delta t_\mathrm{obs}=26.0$-28.7d', 
            transform=ax.transAxes, ha='left', va='top')

    # Formatting
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(7E9, 4E18)
    ax.set_ylim(3E38,6E43)
    ax.legend(loc='lower right')



def flare_sed():
    """ Plot the SED during a flare """

    # Choose the second LRIS flare, since that was simultaneous with 
    # a Chandra observation
    # Add the optical to the plot 
    dat = get_full_opt()
    dat['dt'] = (Time(dat['mjdstart'], format='mjd').jd-vals.t0)
    dt = dat['dt']
    choose_dt = np.logical_and(dt>start_time, dt<end_time)
    # use the NOT gri data
    choose_det = np.logical_and(dat['emag']<99, dat['#instrument']=='NOT/ALFOSC')
    choose = np.logical_and(choose_dt, choose_det)


if __name__=="__main__":
    fig, axarr = plt.subplots(2,1,figsize=(4,5.5))
    ax = axarr[0]
    quiescent_sed(ax)

    ax = axarr[1]
    ax.set_xlabel(
            r"$\nu_\mathrm{obs}$",fontsize=11,
            fontname='sans-serif')

    for ax in axarr:
        ax.tick_params(axis='both', labelsize=11)
        ax.set_ylabel(r"$\nu L_\nu$", fontsize=11,
                fontname='sans-serif')
    plt.tight_layout()
    plt.show()

    #plt.savefig("sed.png", dpi=200, bbox_inches='tight', pad_inches=0.1)
    #plt.close()
