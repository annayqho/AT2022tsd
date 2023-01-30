""" Search for a period in the ULTRASPEC flare data """
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Compute the discrete Fourier transform
def dfft(y, k):
    N = len(y)
    omega = 2*np.pi*k/N
    return sum(y[n]*np.exp(-1j * omega * n) for n in np.arange(N))

# load data
dat = pd.read_fwf(
        "../data/opt/flares_lris_ultraspec.txt", comment='#', 
        names=['MJD','Exp','Filter','Flux','Unc'])

t0 = 59932
dt = dat['MJD'].values-t0
f = dat['Flux'].values
ef = dat['Unc'].values

# Identify the first flare
isdet = f/ef>=5
flare_start = min(dt[isdet])
flare_end = max(dt[np.logical_and(isdet, dt<0.7)])
choose = np.logical_and(dt>=flare_start, dt<=flare_end)
x = dt[choose]
y = f[choose]
ey = ef[choose]

# Get the discrete FT
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

# Identify the second flare
flare_start = min(dt[np.logical_and(isdet, dt>1.72)])
flare_end = max(dt[np.logical_and(isdet, dt<1.79)])
choose = np.logical_and(dt>=flare_start, dt<=flare_end)
x = dt[choose]
y = f[choose]
ey = ef[choose]
plt.errorbar(dt, f, ef, fmt='o', c='lightgrey')
plt.errorbar(dt[isdet], f[isdet], ef[isdet], fmt='o', c='k')
plt.errorbar(x, y, ey, fmt='o', c='red')

# Get the discrete FT
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

#plt.xlim(flare_start-0.01, flare_end+0.01)
