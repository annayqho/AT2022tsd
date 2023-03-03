""" Do time-series analysis of the optical light-curve data """

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time
import sys
import vals
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code/paper_plots")
from get_opt import *
from opt_lc import plot_det, plot_lim


# Generate a power spectrum
def power_spectrum(x,y):
    N = len(x)
    sampling_rate = 1/(x[1]-x[0])
    Y = np.fft.fft(y)/np.sqrt(N)
    f = np.fft.fftfreq(N,1/sampling_rate)
    f,W = f[1:int(N/2)],(2 * np.abs(Y)**2/sampling_rate)[1:int(N/2)]

    # Plot
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.loglog(f,W)
    #ax.set_xlabel("Frequency [1/d]")
    #ax.set_ylabel("power spectral density")
    #plt.show()

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
def get_gband_flare():
    dat = get_full_opt()
    choose = np.logical_and(
            dat['#instrument']=='TNT/ULTRASPEC', dat['flt']=='g')
    dat = dat[choose]
    f = dat['flux'].values
    ef = dat['unc'].values
    dt = dat['mjdstart'].values-dat['mjdstart'].values[0]

    # Parameters
    padding = 0.01
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
def get_rband_flare():
    dat = get_full_opt()
    choose = np.logical_and(
            dat['#instrument']=='TNT/ULTRASPEC', dat['flt']=='r')
    dat = dat[choose]
    f = dat['flux'].values
    ef = dat['unc'].values
    dt = dat['mjdstart'].values-dat['mjdstart'].values[0]

    isdet = f/ef>=5
    flare_start = min(dt[isdet])-0.01
    flare_end = max(dt[isdet])+0.02
    choose = np.logical_and(dt>=flare_start, dt<=flare_end)
    x = dt[choose]
    y = f[choose]
    ey = ef[choose]
    return x,y,ey


def get_rband_noise():
    dat = get_full_opt()
    choose = np.logical_and(
            dat['#instrument']=='TNT/ULTRASPEC', dat['flt']=='r')
    dat = dat[choose]
    f = dat['flux'].values
    ef = dat['unc'].values
    dt = dat['mjdstart'].values-dat['mjdstart'].values[0]
    choose = dt > 0.08
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
    # Construct a periodogram of the ULTRASPEC flares

    # testing things
    #x = np.linspace(0.06, 0.17, 290)
    #y = np.zeros(len(x))
    #y[::5] = np.random.random(len(y[::5]))

    #x,y,ey = get_gband_flare()
    x,y,ey = get_rband_flare()
    sp = np.fft.fft(y)
    freq = np.fft.fftfreq(x.shape[-1])
    dt = x[1]-x[0]
    freq_phys = freq * (1/dt)
    choose = freq_phys > 0
    period_d = 1/freq_phys
    period_m = period_d*24*60
    plt.plot(freq_phys[choose], (sp.real[choose])**2, c='k')

    x,y,ey = get_rband_noise()
    sp = np.fft.fft(y)
    freq = np.fft.fftfreq(x.shape[-1])
    dt = x[1]-x[0]
    freq_phys = freq * (1/dt)
    choose = freq_phys > 0
    period_d = 1/freq_phys
    period_m = period_d*24*60
    plt.plot(freq_phys[choose], (sp.real[choose])**2, c='lightgrey')

    plt.ylabel("|DFFT(flux)|$^2$")
    plt.xlabel("days$^{-1}$")
    plt.tight_layout()

#     # Lomb-Scargle
#     x_input = x_rest * u.second
#     y_input = y * u.uJy
#     ey_input = ey * u.uJy
#     T = max(x_input)-min(x_input)
#     min_freq = 1/T.value
#     max_freq = (1/24)/2
#     freq_step = 1/(10*T.value)
#     freq_grid = np.arange(min_freq, max_freq, step=freq_step) * u.Hz
#     power = LombScargle(x_input, y_input, ey_input).power(freq_grid)
#     plt.plot(1/freq_grid,power)
#     plt.xlabel("Sec")
# 
#     # Plot
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     ax.loglog(pow_x,pow_y)
# 
#     # Get data from the second flare
#     x_obs,y,ey = get_second_flare(dt,f,ef)
#     x_rest = x_obs / (1+vals.z)
#     pow_x,pow_y = power_spectrum(x_rest,y)
# 
#     # Plot
#     ax.loglog(pow_x,pow_y)
#     ax.set_xlabel("Frequency [1/d]")
#     ax.set_ylabel("power spectral density")
#     plt.show()
# # 
# #     # ACF of the noise
# #     x,y,ey = get_noise(dt,f,ef)
# #     acf = autocorr(y)
# #     plt.plot(acf, drawstyle='steps', c='lightgrey')
# # 
#     # ACF of the shorter flare
# 
# 
#     acf = autocorr(y)
#     plt.plot(acf, c='lightblue', drawstyle='steps')
#     plt.errorbar(x, y, ey, fmt='o', c='lightblue')
# # 
# #     # ACF of the longer flare
#     acf = autocorr(y)
#     plt.plot(acf, ls='-', drawstyle='steps', c='darkblue')
# #     plt.tight_layout()
# #     plt.show()
# # 
