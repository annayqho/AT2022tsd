""" Plot showing the radio SED evolution with time """

import numpy as np
import cmasher as cmr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
from scipy.optimize import curve_fit
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
from astropy.time import Time
from get_radio_at2022tsd import *
import vals


def func(x,alpha,y0,x0):
    return y0*(x/x0)**alpha


def get_data():
    """ Get the data for plotting """
    # Load in the data
    dat = get_radio()
    dt = (Time(dat['Date'].values.astype(str),format='isot').jd-vals.t0)/(1+vals.z)
    dat['dt'] = dt

    # Need to add systematic uncertainties to the RMS noise
    # Add 10% to ALMA data
    flux = dat['Flux'].values
    eflux = dat['eFlux'].values
    choose = dat['Tel'].values=='ALMA'
    eflux[choose] = np.sqrt((dat['eFlux'][choose].values)**2+(0.1*flux[choose])**2)

    # Add 5% to Ku band
    choose = np.logical_and(dat['Tel'].values=='VLA', dat['Freq_Obs'].values==15)
    eflux[choose] = np.sqrt((dat['eFlux'][choose].values)**2+(0.05*flux[choose])**2)

    # Add 15% to the three higher bands, for detections
    choose = np.logical_and.reduce((
            dat['Tel'].values=='VLA',dat['Freq_Obs'].values>15,
            dat['Flux']<99))
    eflux[choose] = np.sqrt(
            (dat['eFlux'][choose].values)**2+(0.15*flux[choose])**2)
    
    dat['eFlux'] = eflux
    return dat


def plot_seds(dat, axarr):
    """ Plot a multi-panel SED thing """
    dt = dat['dt']
    flux = dat['Flux']
    eflux = dat['eFlux']

    # Epoch 1
    ax = axarr[0]
    choose = np.logical_and(dt>27, dt<28)
    x = dat['Freq_Obs'][choose]
    y = flux[choose]
    ey = eflux[choose]
    ax.errorbar(x, y, ey, fmt='o', label="$\Delta t=27$-28d", c='k', lw=0.5)
    params, cov = curve_fit(func, x, y, sigma=ey, absolute_sigma=True)
    alpha = params[0]
    y0 = params[1]
    x0 = params[2]
    xvals = np.linspace(10,500)
    yvals = y0*(xvals/x0)**alpha
    ax.plot(
            xvals,yvals,c='lightgrey',ls='-',
            label=r'$f_\nu\propto\nu^{%s}$'%np.round(alpha,1),lw=1)

    # Epoch 2 (ALMA)
    ax= axarr[1]
    choose = np.logical_and(dt>34, dt<37)
    x = dat['Freq_Obs'][choose]
    y = flux[choose]
    ey = eflux[choose]
    ax.errorbar(x, y, ey,
                fmt='o', label="$\Delta t=$34-37d", c='k', lw=0.5)
    params, cov = curve_fit(func, x, y, sigma=ey, absolute_sigma=True)
    alpha = params[0]
    y0 = params[1]
    x0 = params[2]
    xvals = np.linspace(10,500)
    yvals = y0*(xvals/x0)**alpha
    ax.plot(
            xvals,yvals,c='lightgrey',ls='-',
            label=r'$f_\nu\propto\nu^{%s}$'%np.round(alpha,1),lw=1)
     
    # Epoch 3 (NOEMA)
    ax = axarr[2]
    choose = np.logical_and(dt>41, dt<46)
    x = dat['Freq_Obs'][choose]
    y = flux[choose]
    ey = eflux[choose]
    ax.errorbar(x, y, ey, 
                fmt='o', label="$\Delta t=$41-45d", c='k', lw=0.5)
    params, cov = curve_fit(func, x, y, sigma=ey, absolute_sigma=True)
    alpha = params[0]
    y0 = params[1]
    x0 = params[2]
    xvals = np.linspace(10,500)
    yvals = y0*(xvals/x0)**alpha
    ax.plot(
            xvals,yvals,c='lightgrey',ls='-',
            label=r'$f_\nu\propto\nu^{%s}$'%np.round(alpha,1),lw=1)

    choose = np.logical_and(dt>50, dt<51)
    x = dat['Freq_Obs'][choose]
    y = flux[choose]
    ey = eflux[choose]
    ax.errorbar(x, y, ey,
                fmt='D', label="$\Delta t=$50d", c='lightgrey', lw=0.5)
    params, cov = curve_fit(func, x, y, sigma=ey, absolute_sigma=True)
    alpha = params[0]
    y0 = params[1]
    x0 = params[2]
    xvals = np.linspace(10,500)
    yvals = y0*(xvals/x0)**alpha
    ax.plot(
            xvals,yvals,c='lightgrey',ls='--',
            label=r'$f_\nu\propto\nu^{%s}$'%np.round(alpha,1),lw=1)

     
    # Epoch 4 (NOEMA)
    ax = axarr[3]
    choose = np.logical_and(dt>65, dt<70)
    ax.errorbar(dat['Freq_Obs'][choose], flux[choose], eflux[choose], fmt='o',
                label="$\Delta t=$65-70d", c='k')
    #choose = np.logical_and(dt>79, dt<80)
    #ax.errorbar(dat['Freq_Obs'][choose], flux[choose], eflux[choose], fmt='o',
    #            label="$\Delta t=$79d", c='lightgrey')

    # Formatting
    for ax in axarr:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xticks([20,30,50,100,200,500])
        ax.set_yticks([0.02,0.05,0.1, 0.2, 0.5])
        ax.set_xticklabels([20,30,50,100,200,500])
        ax.set_yticklabels([0.02,0.05,0.1, 0.2, 0.5])
        #ax.set_ylabel(r"$f_{\nu}$ (mJy)", fontsize=10)
        ax.set_ylim(0.02,0.7)
        ax.set_xlim(12,500)
        ax.legend(loc='lower right', fontsize=8, handletextpad=0.4,
                  labelspacing=0.1)
    ax = axarr[-1]
    ax.set_xlabel(r"$\nu_\mathrm{obs}$ (GHz)", fontsize=10)


def plot_lc(dat,ax):
    """ Plot the LCs """
    cols = cmr.take_cmap_colors(
            'cmr.freeze', 9, cmap_range=(0.1, 0.9), return_fmt='hex')[::-1]

    nu = dat['Freq_Obs'].astype(int)
    print(np.unique(nu))

    lw = 0.5
    fmt='o'
    ms=3

    for i,val in enumerate([15,22,33,45,77,134,207,350]):
        choose = nu==val

        # Plot the detections
        x = dat['dt'][choose].values
        y = dat['Flux'][choose].values
        ey = dat['eFlux'][choose].values
        isdet = y<99
        ax.errorbar(x[isdet], y[isdet], ey[isdet], 
                    fmt=fmt, c=cols[i], lw=lw, ms=ms)
        ax.plot(x[isdet], y[isdet], lw=1, c=cols[i])

        if val==134:
            ax.text(x[0]/1.02, y[0], str(val),
                    ha='right', va='top',fontsize=8,color=cols[i])
        elif val==15:
            ax.text(x[0]*1.02, y[0], str(val),
                    ha='left', va='top',fontsize=8,color=cols[i])
        elif val==77:
            ax.text(x[0], y[0]*1.02, str(val),
                    ha='center', va='bottom',fontsize=8,color=cols[i])
        else:
            ax.text(x[0]/1.02, y[0], str(val),
                    ha='right', va='center',fontsize=8,color=cols[i])

        # Plot the non-detection
        nondet = y==99
        if sum(nondet)>0:
            lim = ey[nondet][0]*5
            xval = x[nondet][0]
            ax.scatter(xval, lim, marker=fmt, c=cols[i], s=ms*4)
            ax.arrow(xval, lim, 0, -0.015, length_includes_head=True,
                     head_length=0.007, head_width=8, color=cols[i])
            last_x = x[isdet][-1]
            ax.plot(
                    [last_x,xval], [y[isdet][-1],lim], 
                    lw=1, c=cols[i], ls='--')

    ax.axvspan(27,28,color='lightgrey')
    ax.axvspan(34,37,color='lightgrey')
    ax.axvspan(41.2,45.1,color='lightgrey')
    ax.axvspan(65,70,color='lightgrey')
    ax.set_xlabel(r"$\Delta t$ (d)")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xticks([20,30,40,60,80])
    ax.set_yticks([0.02,0.05,0.1, 0.2, 0.5])
    ax.set_xticklabels([20,30,40,60,80])
    ax.set_yticklabels([0.02,0.05,0.1, 0.2, 0.5])
    ax.set_ylabel(r"$f_{\nu}$ (mJy)", fontsize=10)
     

if __name__=="__main__":
    dat = get_data()

    fig = plt.figure(figsize=(6,8))
    gs = gridspec.GridSpec(4,4,hspace=0.2,wspace=0.7)
    ax = fig.add_subplot(gs[1:-1,0:2])
    plot_lc(dat,ax)
    ax1 = fig.add_subplot(gs[0,2:])
    ax2 = fig.add_subplot(gs[1,2:])
    ax3 = fig.add_subplot(gs[2,2:])
    ax4 = fig.add_subplot(gs[3,2:])
    plot_seds(dat, [ax1,ax2,ax3,ax4])

    plt.show()
    #plt.savefig("radio.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    #plt.close()
