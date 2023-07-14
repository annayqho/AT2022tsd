""" Zoom-in of a couple of flares + images for Figure #3 """

import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
from get_opt import *
import vals

# Get the relevant photometry
tab = get_full_opt()

tab_imacs = tab[tab['#instrument']=='Magellan/IMACS']
tab_ultraspec = tab[tab['#instrument']=='TNT/ULTRASPEC']

# Initialize figure
fig,axarr = plt.subplots(2,1,figsize=(7,5))

# Bottom panel: nice plot of the ULTRASPEC r-band flare
ax = axarr[1]
choose = tab_ultraspec['flt']=='r'
toplot = tab_ultraspec[choose]
t0 = toplot['mjdstart'].values[0]
isdet = toplot['sig']>=4
ax.errorbar(
        (toplot['mjdstart'][isdet]-t0)*24*60, toplot['flux'][isdet], 
        toplot['unc'][isdet], fmt='o', c=vals.rc, ms=4, lw=0.5, zorder=5)
ax.errorbar(
        (toplot['mjdstart']-t0)*24*60, toplot['flux'], toplot['unc'],
        fmt='o-',mec=vals.rc,mfc='white',ms=4,lw=0.5, c=vals.rc, zorder=0, 
        alpha=0.4)
ax.axhline(y=0, c='grey', lw=0.5, ls='--')
ax.text(0.01, 0.95, 'ULTRASPEC $r$-band, 30s cadence', ha='left', va='top', 
        transform=ax.transAxes, color=vals.rc)
ax.set_xlim(20, 140)
ax.tick_params(axis='both', labelsize=9)

def new_yaxis(y_i, nueff):
    if y_i > 0:
        y_f = nueff * y_i * 1E-6 * 1E-23 * \
                4 * np.pi * vals.dL_cm**2 / 1E43
    else:
        y_f = 24

ax2 = ax.twinx()
leff = vals.sdss_pivot['r'] # use g-band for all
nueff = 3E18/leff
ymin, ymax = ax.get_ylim()
ax2.set_ylim((new_yaxis(ymin,nueff), new_yaxis(ymax,nueff)))
ax2.plot([],[])
#ax2.set_yscale('log')
ax2.tick_params(axis='both', labelsize=9)
ax2.set_ylabel("$m_\mathrm{AB}$", fontsize=9, rotation=270, labelpad=15)

ax.set_xlabel(r"Minutes after MJD %s" %t0)
ax.set_ylabel(r"Flux Density ($\mu$Jy)")
plt.tight_layout()
plt.show()

