""" Plot the color evolution of the LRIS g+i flare """

import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
from get_opt import *
import vals
from calc_color import *

fig,axarr = plt.subplots(1,3,figsize=(8,2.5))

# Get the data
dat = get_full_opt()
choose = dat['#instrument']=='KeckI/LRIS'
dat = dat[choose]

# Get the g+i flare
choose = np.logical_and(dat['mjdstart']>59871, dat['mjdstart']<59872)
dat = dat[choose]

# Plot flare without Milky Way extinction
ax = axarr[0]
t0 = dat['mjdstart'].values[0]
dt = (dat['mjdstart'].values-t0)*24*60
mag = dat['mag'].values
emag = dat['emag'].values
toplot = dat['flt']=='i'
ax.errorbar(dt[toplot], mag[toplot], emag[toplot], fmt='D', c=vals.ic)
toplot = dat['flt']=='g'
ax.errorbar(dt[toplot], mag[toplot], emag[toplot], fmt='s', c=vals.gc)
ax.invert_yaxis()
ax.text(0.98, 0.98, 'Without MW extinction', transform=ax.transAxes, ha='right', va='top', fontsize=8)
ax.set_ylabel("Apparent Mag")
ax.set_xlabel("Minutes")

# With Milky Way extinction
ax = axarr[1]
toplot = dat['flt']=='i'
ax.errorbar(dt[toplot], mag[toplot]-vals.ext['i'], emag[toplot], fmt='D', c=vals.ic)
toplot = dat['flt']=='g'
ax.errorbar(dt[toplot], mag[toplot]-vals.ext['g'], emag[toplot], fmt='s', c=vals.gc)
ax.invert_yaxis()
ax.text(0.98, 0.98, 'With MW extinction', transform=ax.transAxes, ha='right', va='top', fontsize=8)
ax.set_xlabel("Minutes")

# Spectral index
ax = axarr[2]
beta, ebeta = calc_color_spindex(dat[choose], 'g', 'i')
ax.errorbar(dt[toplot], beta, ebeta, c='grey', fmt='o')
ax.set_ylabel(r"$\beta$ where $f_\nu \propto \nu^{\beta}$")
ax.set_xlabel("Minutes")

plt.tight_layout()
plt.savefig("lris_flare_color_evolution.png", dpi=200)
plt.close()
#plt.show()

