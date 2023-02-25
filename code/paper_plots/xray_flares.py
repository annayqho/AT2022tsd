""" Plot the Chandra X-ray flares """

import numpy as np
import pandas as pd
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt

dd = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/xray"

# Read the Chandra summary
dat = pd.read_csv(dd + "/chandra_flux_summary.txt")
oids = dat['OBSID'].values.astype(str)
dates = dat['Start'].values

# Sort in order of dates
order = np.argsort(dates)
oids = oids[order]
dates = dates[order]

# Initialize the figure...
fig,axarr = plt.subplots(4,2,figsize=(8,6))

for i,oid in enumerate(oids):
    ax = axarr.flatten()[i]
    ff = dd + "/" + oid + "/repro/xray_flare_lc.txt"
    bf = dd + "/" + oid + "/repro/xray_bkg_lc.txt"
    ff2 = dd + "/" + oid + "/repro/xray_src_offset_lc.txt"

    dat = np.loadtxt(ff)
    ax.errorbar(
            dat[:,0]/60, dat[:,1], yerr=dat[:,2], marker="o", c='k',
            label='AT2022tsd',zorder=2,ms=4, lw=0.5)

    # Scale the count rate by the ratio of the size of the background
    # to the size of the source
    src_size = 2
    with open(dd + "/" + oid + "/repro/bkg_offset.reg") as inputf:
        l = inputf.readlines()[-1]
        bkg_size = float(l.split(',')[-1][0:-2])/2

    dat = np.loadtxt(bf)
    ax.plot(dat[:,0]/60, dat[:,1]*(src_size**2/bkg_size**2), c='lightgrey', 
            lw=2, label='Background', zorder=0)

    if i==0:
        ax.legend(loc='upper left', fontsize=8)

    dat = dates[i].split('T')[0]
    if np.logical_or(i==0, i==2):
        ax.text(
                0.88, 0.95, dat, ha='right', va='top', 
                fontsize=8, transform=ax.transAxes)
    else:
        ax.text(
                0.98, 0.95, dat, ha='right', va='top', 
                fontsize=8, transform=ax.transAxes)

axarr.flatten()[-1].axis('off')

for ax in axarr[:,0]:
    ax.set_ylabel("Count Rate")
for ax in axarr[-1,:]:
    ax.set_xlabel("Minutes")

plt.tight_layout()
plt.savefig("xray_flare_lc.png", dpi=200)
plt.close()
#plt.show()
