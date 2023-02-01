""" Plot the ULTRASPEC flares """

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import matplotlib.patches as patches
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["pdf.fonttype"] = 42
from astropy.time import Time
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
from get_opt import get_ipac,get_flares,get_ultraspec
import vals

def plot_ultraspec_panel(ax, dat, filt, m, col, zoom=False, inset=False):
    """ Plot a light-curve panel """
    choose = dat['Filt']==filt
    x = dat['MJD'][choose]
    y = dat['Flux'][choose]
    ey = dat['Unc'][choose]
    t0 = min(x)

    s=2
    l=0.2

    # Definition of a detection: 5-sigma
    # Otherwise you get what seem to be spurious points
    det = y/ey >= 5

    if zoom:
        t0 = min(x[det]*24*60)-5
        dt = x*24*60-t0
    else:
        dt = (x-t0)*24

    # Plot detections
    ax.errorbar(dt[det], y[det], ey[det],
                c=col, fmt=m, ms=s,lw=l)

    # Plot the non-detections
    ax.errorbar(dt[~det], y[~det], ey[~det],lw=l,
                mec=col, fmt=m, ms=s, alpha=0.3, mfc='white', c=col)

    # Plot a y=0 line
    ax.axhline(y=0, c='grey', lw=0.5)

    if zoom:
        if filt=='g':
            ax.set_xlim(35, max(dt[det])+5)
        else:
            ax.set_xlim(min(dt[det])-5, max(dt[det])+5)

    if inset:
        ax.set_xlim(0, 10)

    if inset==False and zoom==False:
        ax.set_xlim(min(dt), max(dt))

    ax.set_ylim(min(y)-1, max(y)+3)


# Initialize figure
fig,axarr = plt.subplots(3,2,figsize=(7,6))

# Plot the Magellan flare
tel,mjd,filt,mag,emag,limmag,flare = get_flares()
choose = tel=='Magellan'
ax = axarr[0,0]
y = (10**((mag[choose]-8.9)/(-2.5))) * 1E6
ey = y*emag[choose]*(np.log(10))/2.5
ax.errorbar((mjd[choose]-mjd[choose][0])*24*60, y, ey,
            fmt='s', c=vals.gc)
ax.text(0.02, 0.98, 'IMACS $g$-band', transform=ax.transAxes,
        ha='left', va='top', fontsize=8)
ax.text(0.02, 0.88, '2022-12-15', transform=ax.transAxes,
        ha='left', va='top', fontsize=8)

# Plot the LT flare
choose = np.logical_and(tel=='LT', flare=='*')
ax = axarr[0,1]
y = (10**((mag[choose]-8.9)/(-2.5))) * 1E6
ey = y*emag[choose]*(np.log(10))/2.5
ax.errorbar((mjd[choose]-mjd[choose][0])*24*60, y, ey,
            fmt='s', c=vals.gc)
ax.text(0.98, 0.98, 'LT $g$-band', transform=ax.transAxes,
        ha='right', va='top', fontsize=8)
ax.text(0.98, 0.88, '2022-12-16', transform=ax.transAxes,
        ha='right', va='top', fontsize=8)

# Get the ULTRASPEC data
dat = get_ultraspec()

# Top-left panel: r-band flare 
ax = axarr[1,0]
plot_ultraspec_panel(ax, dat, 'r', 'o', vals.rc)
ax.text(0.98, 0.95, 'ULTRASPEC $r$-band', transform=ax.transAxes,
        ha='right', va='top', fontsize=8)
ax.text(0.98, 0.85, '2022-12-19', transform=ax.transAxes,
        ha='right', va='top', fontsize=8)

# Top-right panel: zoom in
plot_ultraspec_panel(axarr[1,1], dat, 'r', 'o', vals.rc, zoom=True)

# Bottom-left panel: g-band flare
ax = axarr[2,0]
plot_ultraspec_panel(ax, dat, 'g', 's', vals.gc)
ax.text(0.98, 0.1, 'ULTRASPEC $g$-band', transform=ax.transAxes,
        ha='right', va='top', fontsize=8)
ax.text(0.98, 0.17, '2022-12-20', transform=ax.transAxes,
        ha='right', va='top', fontsize=8)

# Zoom-in of the first flare
axins = axarr[2,0].inset_axes([0.03, 0.55, 0.35, 0.4])
axins.set_yticks([])
plot_ultraspec_panel(axins, dat, 'g', 's', vals.gc, inset=True, zoom=True)
axins.tick_params(axis='x', labelsize=8, pad=0.5)
axins.set_xlabel("Minutes", fontsize=8, labelpad=1)
axins.set_ylim(-1, 15)

# Bottom-right panel: zoom in of the second flare
plot_ultraspec_panel(axarr[2,1], dat, 'g', 's', vals.gc, zoom=True)

# Formatting
axarr[2,0].set_xlabel("Hours")
axarr[2,1].set_xlabel("Minutes")
for ax in axarr[:,0]:
    ax.set_ylabel("Flux Density ($\mu$Jy)")

#plt.savefig(
#        "ultraspec_flares.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
#plt.close()
plt.show()


