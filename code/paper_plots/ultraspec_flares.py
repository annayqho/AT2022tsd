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

def plot_panel(ax, dat, filt, m, col, zoom=False, inset=False):
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
fig,axarr = plt.subplots(2,2,figsize=(7,5))

# Get the data
dat = get_ultraspec()

# Top-left panel: r-band flare 
plot_panel(axarr[0,0], dat, 'r', 'o', vals.rc)

# Top-right panel: zoom in
plot_panel(axarr[0,1], dat, 'r', 'o', vals.rc, zoom=True)

# Bottom-left panel: g-band flare
plot_panel(axarr[1,0], dat, 'g', 's', vals.gc)

# Zoom-in of the first flare
axins = axarr[1,0].inset_axes([0.03, 0.55, 0.35, 0.4])
axins.set_yticks([])
plot_panel(axins, dat, 'g', 's', vals.gc, inset=True, zoom=True)
axins.tick_params(axis='x', labelsize=8, pad=0.5)
axins.set_xlabel("Minutes", fontsize=8, labelpad=1)
axins.set_ylim(-1, 15)

# Bottom-right panel: zoom in of the second flare
plot_panel(axarr[1,1], dat, 'g', 's', vals.gc, zoom=True)

# Formatting
axarr[1,0].set_xlabel("Hours")
axarr[1,1].set_xlabel("Minutes")
for ax in axarr[:,0]:
    ax.set_ylabel("Flux Density ($\mu$Jy)")

#plt.savefig(
#        "ultraspec_flares.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
#plt.close()
plt.show()


