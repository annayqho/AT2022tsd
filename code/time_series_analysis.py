""" Search for a period in the ULTRASPEC flare data """
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
    isdet = f/ef>=5
    flare_start = min(dt[isdet])-0.01
    flare_end = max(dt[np.logical_and(isdet, dt<0.7)])+0.01
    choose = np.logical_and(dt>=flare_start, dt<=flare_end)
    x = dt[choose]
    y = f[choose]
    ey = ef[choose]
    return x,y,ey


# Get the second flare
def get_second_flare(dt,f,ef):
    isdet = f/ef>=5
    flare_start = min(dt[np.logical_and(isdet, dt>1.72)])-0.01
    flare_end = max(dt[np.logical_and(isdet, dt<1.79)])+0.01
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
def period_search(y):
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


if __name__=="__main__":
    # load data
    dat = pd.read_fwf(
            "../data/opt/flares_lris_ultraspec.txt", comment='#',
            names=['MJD','Exp','Filter','Flux','Unc'])
    t0 = 59932
    dt = dat['MJD'].values-t0
    f = dat['Flux'].values
    ef = dat['Unc'].values
    #plt.errorbar(dt, f, ef, fmt='o', c='lightgrey')

    # ACF of the noise
    x,y,ey = get_noise(dt,f,ef)
    acf = autocorr(y)
    plt.plot(acf, drawstyle='steps', c='lightgrey')

    # ACF of the shorter flare
    x,y,ey = get_first_flare(dt,f,ef)
    acf = autocorr(y)
    plt.plot(acf, c='lightblue', drawstyle='steps')
    #plt.errorbar(x, y, ey, fmt='o', c='lightblue')

    # ACF of the longer flare
    x,y,ey = get_second_flare(dt,f,ef)
    acf = autocorr(y)
    plt.plot(acf, ls='-', drawstyle='steps', c='darkblue')
    plt.tight_layout()
    plt.show()

