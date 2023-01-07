""" Plot the VLA + NOEMA + SMA data """

import matplotlib.pyplot as plt
import pandas as pd
from nufnu_time import *

# two panels
fig,axarr = plt.subplots(1,2,figsize=(7,4))

# Left panel: the mm light curve in context
ax = axarr[0]
run(ax)
ax.set_yscale('log')
ax.set_xscale('log')

# Right panel: the radio/mm SED
ax = axarr[1]

dat = pd.read_csv("../data/radio.txt", delimiter=',')
dt = dat['dt'].values
nu = dat['Freq_Obs'].values
f = dat['Flux_mJy'].values
ef = dat['eFlux'].values

choose = np.logical_and(np.logical_and(dt>24,dt<28), dat['Tel']!='SMA')
ax.errorbar(
        nu[choose], f[choose].astype(float), ef[choose].astype(float), 
        fmt='o', c='grey', label=r'$\Delta t=$26d')
choose = np.logical_and(dt>32,dt<36)
ax.errorbar(
        nu[choose], f[choose].astype(float), ef[choose].astype(float), 
        fmt='s', c='k', label=r'$\Delta t=$34d')
choose = np.logical_and(dt>41,dt<46)
ax.errorbar(
        nu[choose], f[choose].astype(float), ef[choose].astype(float), 
        fmt='D', c='purple', label=r'$\Delta t=$43d')
choose = np.logical_and(dt>50,dt<54)
ax.errorbar(
        nu[choose], f[choose].astype(float), ef[choose].astype(float), 
        fmt='X', c='green', label=r'$\Delta t=$51d')

xvals = np.linspace(5,400)
yvals = 0.1*(xvals/40)**2
#ax.plot(
#        xvals, yvals, c='k', ls='--', lw=0.5, 
#        label=r'$f_\nu \propto \nu^2$')
#ax.set_title(r"$\Delta t_\mathrm{obs}\approx$7 d")
ax.text(20,1E-1,'VLA',color='k')
ax.text(100,0.6,'NOEMA',color='k',ha='left',va='top')
#ax.axvspan(211,275,color='darkgrey')
#ax.axvspan(275,373,color='lightgrey')
#ax.axvspan(385,500,color='darkgrey')
#ax.axvspan(602,720,color='lightgrey')
ax.legend(loc='upper left')

ax.set_ylim(0.02, 0.8)
ax.set_xlim(8,450)
for ax in axarr:
    ax.set_xlabel(r"$\nu_\mathrm{obs}$ (GHz)")
axarr[0].set_ylabel(
            r"$L_{\nu}$ (erg$\,$s$^{-1}$Hz$^{-1}$)")
ax.set_ylabel(r"$F_\nu$ (mJy)")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xticks([20,30,50,100,200,400])
ax.set_xticklabels([20,30,50,100,200,400])
ax.set_yticks([0.03, 0.05, 0.1, 0.2, 0.5])
ax.set_yticklabels([0.03, 0.05, 0.1, 0.2, 0.5])


# final formatting
plt.tight_layout()
plt.show()
#plt.savefig("mm_radio.png", dpi=200)
#plt.close()
