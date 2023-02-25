""" Plot the Chandra X-ray flares """

import numpy as np
import pandas as pd
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt

dd = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/xray"
oids = ['26641', '26642', '26643', '26644', '26645', '27639', '27643']
oids = ['26641', '26642']

# Initialize the figure...
fig,axarr = plt.subplots(3,2,figsize=(8,6))

for i,oid in enumerate(oids):
    ax = axarr.flatten()[i]
    ff = dd + "/" + oid + "/repro/xray_flare_lc.txt"
    bf = dd + "/" + oid + "/repro/xray_bkg_lc.txt"
    ff2 = dd + "/" + oid + "/repro/xray_src_offset_lc.txt"

    dat = np.loadtxt(ff)
    ax.errorbar(
            dat[:,0], dat[:,1], yerr=dat[:,2], marker="o", c='k',
            label='AT2022tsd, background subtracted',zorder=2,ms=5)

    #dat = np.loadtxt(ff2)
    #ax.errorbar(
    #        dat[:,0], dat[:,1], yerr=dat[:,2], marker="o", c='grey',
    #        label='2\'\'-radius nearby region',zorder=0)

    dat = np.loadtxt(bf)
    ax.plot(dat[:,0], dat[:,1]*(2**2/10**2), c='lightgrey', lw=2, 
            label='Background, scaled to source-region area', zorder=0)
    #ax.errorbar(
    #        dat[:,0], dat[:,1], yerr=dat[:,2], marker="o", c='lightgrey',
    #        label='Background',zorder=0)

    if i==0:
        ax.legend(loc='upper left', fontsize=8)

for ax in axarr[:,0]:
    ax.set_ylabel("Count Rate")
for ax in axarr[-1,:]:
    ax.set_xlabel("Seconds")

plt.tight_layout()
#plt.savefig("first_flare_lc.png", dpi=200)
#plt.close()
plt.show()
