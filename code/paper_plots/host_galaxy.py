""" Figure showing the position of the transient in its host galaxy,
as well as the position of the host galaxy in the M*-SFR plane """

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import make_lupton_rgb

fig,axarr = plt.subplots(1,2,figsize=(8,4))

# Position of the transient in the host galaxy
ax = axarr[0]
ddir = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/host/"
g = fits.open(
        ddir + "cutout_rings.v3.skycell.1513.011.stk.g.unconv.fits")[0].data
r = fits.open(
        ddir + "cutout_rings.v3.skycell.1513.011.stk.r.unconv.fits")[0].data
i = fits.open(
        ddir + "cutout_rings.v3.skycell.1513.011.stk.i.unconv.fits")[0].data
z = fits.open(
        ddir + "cutout_rings.v3.skycell.1513.011.stk.z.unconv.fits")[0].data
rgb = make_lupton_rgb(i/4, r/2, g, Q=10, stretch=0.5) #, minimum=[2,2,2], stretch=5, Q=8)
ax.imshow(rgb[35:240,35:240,:], origin='lower')

#ax = axarr[1]

plt.show()
