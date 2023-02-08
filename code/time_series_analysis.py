""" Do time-series analysis of the optical light-curve data """

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time
import sys
import vals
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code/paper_plots")
from get_opt import get_ipac,get_flares
from opt_lc import plot_det, plot_lim


# Generate a power spectrum
def power_spectrum(x,y):
    N = len(x)
    sampling_rate = 1/(x[1]-x[0])
    Y = np.fft.fft(y)/np.sqrt(N)
    f = np.fft.fftfreq(N,1/sampling_rate)
    f,W = f[1:int(N/2)],(2 * np.abs(Y)**2/sampling_rate)[1:int(N/2)]

    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.loglog(f,W)
    ax.set_xlabel("Frequency [1/d]")
    ax.set_ylabel("power spectral density")
    plt.show()

    return f,W


# Compute the discrete Fourier transform
def dfft(y, k):
    N = len(y)
    omega = 2*np.pi*k/N
    return sum(y[n]*np.exp(-1j * omega * n) for n in np.arange(N))


#  Compute the autocorrelation
def autocorr(y):
    yfft = np.fft.fft(y)
    acf = np.fft.fftshift(np.fft.ifft(yfft * np.conj(yfft)))
    return acf


# Get the first flare
def get_first_flare(dt,f,ef):
    # Parameters
    padding = 0.01
    padding = 0
    sig = 4

    isdet = f/ef>=sig
    flare_start = min(dt[isdet])-padding
    flare_end = max(dt[np.logical_and(isdet, dt<0.7)])+padding
    choose = np.logical_and(dt>=flare_start, dt<=flare_end)
    x = dt[choose]
    y = f[choose]
    ey = ef[choose]
    return x,y,ey


# Get the second flare
def get_second_flare(dt,f,ef):
    isdet = f/ef>=5
    flare_start = min(dt[np.logical_and(isdet, dt>1.72)])-0.01
    flare_end = max(dt[np.logical_and(isdet, dt<1.79)])+0.02
    choose = np.logical_and(dt>=flare_start, dt<=flare_end)
    x = dt[choose]
    y = f[choose]
    ey = ef[choose]
    return x,y,ey


def get_noise(dt,f,ef):
    choose = np.logical_and(dt>=0.7, dt<=0.8)
    x = dt[choose]
    y = f[choose]
    ey = ef[choose]
    return x,y,ey


# Get the discrete FT
def period_search_dfft(y):
    """ Search for a period within a single flare """
    Nbin = len(y)/2
    k = np.arange(Nbin)
    Y = np.array([np.abs(dfft(y, kval)) for kval in k])
    bins = np.array([kval*0.5 for kval in k])
    plt.bar(bins, Y, width=0.5)
    plt.axvline(x=1, ls='--', c='k')
    plt.xlabel("Minutes")
    plt.ylabel("Amplitude")
    plt.tight_layout()
    plt.show()


def period_search():
    """ Search for a periodicity from one flare to the next """
    # Impose a threshold magnitude
    thresh = 21 # flare needs to be brighter than 21 mag

    fig,ax = plt.subplots(1,1,figsize=(6,3.5))

    # Plot the main transient LC
    jd,filt,mag,emag,fujy,efujy = get_ipac()
    dt = jd-vals.t0

    # Just get the main LC, not the flares
    choose = dt < 20
    dt = dt[choose]
    jd = jd[choose]
    filt = filt[choose]
    mag = mag[choose]
    emag = emag[choose]

    # Plot the g-band detections
    choose = np.logical_and(filt=='g', emag<99)
    plot_det(
            ax, dt[choose], mag[choose], emag[choose],
            'g', lines=True)

    # Plot the g-band limits
    choose = np.logical_and(filt=='g', emag==99)
    plot_lim(ax, dt[choose], mag[choose], 'g')

    # Plot the r-band detections
    choose = np.logical_and(filt=='r', emag<99)
    plot_det(ax, dt[choose], mag[choose], emag[choose], 'r',
             lines=True)

    # Plot the r-band limits
    choose = np.logical_and(filt=='r', emag==99)
    plot_lim(ax, dt[choose], mag[choose], 'r')

    # Now plot the flares
    # First, ZTF flares
    jd,filt,mag,emag,fujy,efujy = get_ipac()
    dt = jd-vals.t0

    # Just get the flares
    choose = dt > 20
    dt = dt[choose]
    jd = jd[choose]
    filt = filt[choose]
    mag = mag[choose]
    emag = emag[choose]

    # For plotting purposes
    filt[filt=='sdssg'] = ['g']*sum(filt=='sdssg')
    filt[filt=='sdssr'] = ['r']*sum(filt=='sdssr')

    # Plot the g-band detections
    choose = np.logical_and(filt=='g', emag<99)
    plot_det(
            ax, dt[choose], mag[choose], emag[choose], 'g')

    # Plot the g-band limits
    choose = np.logical_and(filt=='g', emag==99)
    plot_lim(ax, dt[choose], mag[choose], 'g')

    # Plot the r-band detections
    choose = np.logical_and(filt=='r', emag<99)
    plot_det(ax, dt[choose], mag[choose], emag[choose], 'r')

    # Plot the r-band limits
    choose = np.logical_and(filt=='r', emag==99)
    plot_lim(ax, dt[choose], mag[choose], 'r')

    # Non-ZTF flares
    tel,mjd,filt,mag,emag,limmag,flare = get_flares()
    jd = Time(mjd, format='mjd').jd
    dt = jd-vals.t0

    # For plotting purposes
    filt[filt=='sdssg'] = ['g']*sum(filt=='sdssg')
    filt[filt=='sdssr'] = ['r']*sum(filt=='sdssr')

    # Plot the g-band flares
    choose = np.logical_and(flare=='*', filt=='g')
    plot_det(ax,dt[choose],mag[choose],emag[choose],'g')

    # Plot the g-band limits
    choose = np.logical_and(filt=='g', emag==99)
    plot_lim(ax, dt[choose], mag[choose], 'g')

    # Plot the r-band flares
    choose = np.logical_and(flare=='*', filt=='r')
    plot_det(ax,dt[choose],mag[choose],emag[choose],'r')

    # Plot the r-band limits
    choose = np.logical_and(filt=='r', emag==99)
    plot_lim(ax, dt[choose], mag[choose], 'r')

    plt.gca().invert_yaxis()
    plt.ylim(22, 20)
    plt.axhline(y=21, c='k')
    plt.show()



if __name__=="__main__":
    #period_search()
    # load data
    dat = pd.read_fwf(
            "../data/opt/flares_lris_ultraspec.txt", comment='#',
            names=['MJD','Exp','Filter','Flux','Unc'])
    t0 = 59932
    dt = dat['MJD'].values-t0
    f = dat['Flux'].values
    ef = dat['Unc'].values
    plt.errorbar(dt, f, ef, fmt='o', c='lightgrey')
# 
#     # ACF of the noise
#     x,y,ey = get_noise(dt,f,ef)
#     acf = autocorr(y)
#     plt.plot(acf, drawstyle='steps', c='lightgrey')
# 
    # ACF of the shorter flare
    x,y,ey = get_first_flare(dt,f,ef)
    acf = autocorr(y)
    plt.plot(acf, c='lightblue', drawstyle='steps')
    plt.errorbar(x, y, ey, fmt='o', c='lightblue')
# 
#     # ACF of the longer flare
    x,y,ey = get_second_flare(dt,f,ef)
    acf = autocorr(y)
    plt.plot(acf, ls='-', drawstyle='steps', c='darkblue')
#     plt.tight_layout()
#     plt.show()
# 
