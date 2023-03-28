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


def plot_seds(dat, ax):
    """ Plot a multi-panel SED thing """
    cols = cmr.take_cmap_colors(
            'cmr.ember', 6, cmap_range=(0.1, 0.9), return_fmt='hex')[::-1]
    
    dt = dat['dt']
    flux = dat['Flux']
    eflux = dat['eFlux']

    mins = [26,33,40.2,64,112,141]
    maxs = [27,36,44.2,69,113,157]
    markers = ['o', 's', 'D', '>', '*', '<']
    sizes = [5, 5, 5, 6, 10, 6]

    for i,minval in enumerate(mins):
        choose = np.logical_and(dt>minval, dt<maxs[i])
        x = dat['Freq_Obs'][choose].values*(1+vals.z)
        order = np.argsort(x)
        x = x[order]
        y = flux[choose].values[order]
        ey = eflux[choose].values[order]

        if minval==141:
            print(x, y, ey)

        # Plot detections
        isdet = y<99
        if minval==141:
            print(x[isdet], y[isdet], ey[isdet])
        ax.errorbar(
                x[isdet], y[isdet], ey[isdet], 
                fmt=markers[i], c=cols[i], lw=0.5, ms=sizes[i],
                label="$\Delta t=%s$-%sd" %(minval,maxs[i]))
        ax.plot(x[isdet],y[isdet],c=cols[i],lw=2)

        # Plot the nondetection
        if sum(~isdet)>0:
            print("There are non-detections")
            xval = x[~isdet][0]
            lims = ey[~isdet]*3
            yval = lims[0] # 3-sigma limits
            ax.scatter(x[~isdet], lims, marker=markers[i], 
                       edgecolor=cols[i], facecolor='white',
                       s=sizes[i]*10, zorder=100)
            #ax.arrow(x[~isdet][0],ey[~isdet][0]*5,0,-0.015,
            #         length_includes_head=True,
            #         head_length=0.007, head_width=8, color=cols[i],zorder=10)
            ax.plot([x[isdet][-1],xval], [y[isdet][-1],yval], 
                    c=cols[i], lw=2, ls='--', zorder=10)
            ax.plot([x[isdet][-1],x[~isdet][-1]], [y[isdet][-1],lims[-1]], 
                    c=cols[i], lw=2, ls='--', zorder=10)
            #ax.plot(x[~isdet], lims, c=cols[i], lw=2, ls='--', zorder=10)

    #params, cov = curve_fit(func, x, y, sigma=ey, absolute_sigma=True)
    #alpha = params[0]
    #y0 = params[1]
    #x0 = params[2]
    xvals = np.linspace(10,500)
    alpha = 5/2
    x0 = 50
    y0 = 0.3
    yvals = y0*(xvals/x0)**alpha
    ax.plot(xvals,yvals,lw=0.5,c='grey')
    ax.text(60,0.5,r'$f_\nu\propto\nu^{5/2}$',c='grey',ha='right',fontsize=8)

    alpha = 1
    x0 = 30
    y0 = 0.06
    yvals = y0*(xvals/x0)**alpha
    ax.plot(xvals,yvals,lw=1,c='grey',ls=':')
    ax.text(300,0.5,r'$f_\nu\propto\nu^{1}$',c='grey',ha='left',fontsize=8)

    # Formatting
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xticks([20,30,50,100,200,500])
    ax.set_yticks([0.02,0.05,0.1, 0.2, 0.5])
    ax.set_xticklabels([20,30,50,100,200,500])
    ax.set_yticklabels([0.02,0.05,0.1, 0.2, 0.5])
    #ax.set_ylabel(r"$f_{\nu}$ (mJy)", fontsize=10)
    ax.set_ylim(0.02,0.7)
    ax.set_xlim(1,600)
    ax.legend(loc='lower right', fontsize=8, handletextpad=0.4,
              labelspacing=0.1)
    ax.set_xlabel(r"$\nu_\mathrm{rest}$ (GHz)", fontsize=10)


def plot_lc(dat,ax):
    """ Plot the LCs """
    cols = cmr.take_cmap_colors(
            'cmr.cosmic', 9, cmap_range=(0.1, 0.9), return_fmt='hex')[::-1]

    nu = dat['Freq_Obs'].astype(int)

    lw = 0.5
    fmt='o'
    ms=3

    for i,val in enumerate([15,22,33,45,77,134,207,350]):
        choose = nu==val

        # Plot the detections
        x = dat['dt'][choose].values
        y = dat['Flux'][choose].values
        ey = dat['eFlux'][choose].values
        isdet = np.logical_and(y<99, y/ey>3)
        ax.errorbar(x[isdet], y[isdet], ey[isdet], 
                    fmt=fmt, c=cols[i], lw=lw, ms=ms, zorder=10)
        ax.plot(x[isdet], y[isdet], lw=1, c=cols[i], zorder=10)

        if val==134:
            ax.text(x[0]/1.02, y[0], str(val) + 'GHz',
                    ha='right', va='top',fontsize=8,color=cols[i])
        elif val==15:
            ax.text(x[0]*1.02, y[0], str(val) + 'GHz',
                    ha='left', va='top',fontsize=8,color=cols[i])
        elif val==77:
            ax.text(x[-1], ey[-1]*3*1.02, str(val) + 'GHz',
                    ha='right', va='bottom',fontsize=8,color=cols[i])
        else:
            ax.text(x[0]/1.02, y[0], str(val) + 'GHz',
                    ha='right', va='center',fontsize=8,color=cols[i])

        # Plot the non-detection
        nondet = ~isdet
        if sum(nondet)>0:
            lim = ey[nondet][0]*3
            xval = x[nondet][0]
            ax.scatter(
                    xval, lim, marker=fmt, edgecolor=cols[i], 
                    facecolor='white', s=ms*6, zorder=20)
            #ax.arrow(xval, lim, 0, -lim/5, length_includes_head=True,
            #         head_length=lim/10, head_width=val/10, color=cols[i], zorder=10)
            last_x = x[isdet][-1]
            ax.plot(
                    [last_x,xval], [y[isdet][-1],lim], 
                    lw=1, c=cols[i], ls='--', zorder=10)

    ax.axvspan(27,28,color='lightgrey')
    ax.axvspan(34,37,color='lightgrey')
    ax.axvspan(41.2,45.1,color='lightgrey')
    ax.axvspan(65,70,color='lightgrey')
    ax.axvspan(113,114,color='lightgrey',zorder=0)
    ax.set_xlabel(r"$\Delta t_\mathrm{rest}$ (d)")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(18,180)
    ax.set_xticks([20,30,40,60,80,120,180])
    ax.set_xticklabels([20,30,40,60,80,120,180])
    ax.set_yticks([0.02,0.05,0.1, 0.2, 0.5])
    ax.set_yticklabels([0.02,0.05,0.1, 0.2, 0.5])
    ax.set_ylabel(r"$f_{\nu}$ (mJy)", fontsize=10)
     

if __name__=="__main__":
    dat = get_data()

    fig,axarr = plt.subplots(1,2, figsize=(8,4))
    ax = axarr[0]
    plot_lc(dat,ax)
    ax = axarr[1]
    plot_seds(dat,ax)

    plt.show()
    #plt.savefig("radio.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    #plt.close()
