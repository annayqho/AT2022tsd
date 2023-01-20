""" Plot showing the radio SED evolution with time """

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
from astropy.time import Time
from get_radio_at2022tsd import *
import vals

fig,axarr = plt.subplots(4,1,figsize=(6,8))

# Load in the data
dat = get_radio()
dt = (Time(dat['Date'].values.astype(str),format='isot').jd-vals.t0)/(1+vals.z)
dat['dt'] = dt
print(dat)

# Need to add systematic uncertainties to the RMS noise
# Add 10% to ALMA data
flux = dat['Flux'].values
eflux = dat['eFlux'].values
choose = dat['Tel'].values=='ALMA'
eflux[choose] = np.sqrt((dat['eFlux'][choose].values)**2+(0.1*flux[choose])**2)

# Epoch 1
ax = axarr[0]
choose = np.logical_and(dt>27, dt<28)
x = dat['Freq_Obs'][choose]
y = flux[choose]
ey = eflux[choose]
ax.errorbar(x, y, ey, fmt='o', label="$\Delta t=27$-28d", c='k', lw=0.5)
xvals = np.linspace(10,300)
yvals = (0.1)*(xvals/40)**1
plt.plot(xvals,yvals,c='Crimson',ls='-',label=r'$f_\nu\propto\nu^{1}$',lw=0.5)
yvals = (0.1)*(xvals/40)**2.5
plt.plot(xvals,yvals,c='Crimson',ls='--',label=r'$f_\nu\propto\nu^{5/2}$',lw=0.5)
yvals = (0.3)*(xvals/200)**(-1)
plt.plot(xvals,yvals,c='Crimson',ls=':',label=r'$f_\nu\propto\nu^{-1}$',lw=0.5)

# Epoch 2 (ALMA)
ax= axarr[1]
choose = np.logical_and(dt>34, dt<37)
ax.errorbar(dat['Freq_Obs'][choose], flux[choose], eflux[choose], fmt='s',
            label="$\Delta t=$34-37d", c='darkgrey', lw=0.5)
# 
# Epoch 3 (NOEMA)
ax = axarr[2]
choose = np.logical_and(dt>41, dt<46)
ax.errorbar(dat['Freq_Obs'][choose], flux[choose], eflux[choose], fmt='D',
            label="$\Delta t=$41.2-45.1d", c='k', lw=0.5)
choose = np.logical_and(dt>50, dt<51)
ax.errorbar(dat['Freq_Obs'][choose], flux[choose], eflux[choose], fmt='D',
            label="$\Delta t=$50d", c='lightgrey', lw=0.5)

 
# Epoch 4 (NOEMA)
ax = axarr[3]
choose = np.logical_and(dt>65, dt<70)
ax.errorbar(dat['Freq_Obs'][choose], flux[choose], eflux[choose], fmt='o',
            label="$\Delta t=$65-70d", c='lightgrey')
choose = np.logical_and(dt>79, dt<80)
ax.errorbar(dat['Freq_Obs'][choose], flux[choose], eflux[choose], fmt='o',
            label="$\Delta t=$79d", c='k')

# Formatting
for ax in axarr:
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xticks([20,30,50,100,200,500])
    ax.set_yticks([0.03,0.05,0.1, 0.2, 0.5])
    ax.set_xticklabels([20,30,50,100,200,500])
    ax.set_yticklabels([0.03,0.05,0.1, 0.2, 0.5])
    ax.set_ylabel(r"$f_{\nu}$ (mJy)", fontsize=10)
    ax.set_ylim(0.02,0.7)
    ax.set_xlim(12,500)
    ax.legend(loc='upper left', fontsize=8)
ax = axarr[-1]
ax.set_xlabel(r"$\nu_\mathrm{obs}$ (GHz)", fontsize=10)

plt.tight_layout()
plt.show()

 
