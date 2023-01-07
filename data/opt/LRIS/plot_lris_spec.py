import matplotlib.pyplot as plt
import sys
import numpy as np
import pandas as pd
import astropy.units as u
from specutils import Spectrum1D
from specutils.manipulation import (box_smooth, gaussian_smooth, trapezoid_smooth)
from astropy.nddata import StdDevUncertainty

dat = pd.read_fwf("lris_20221006.ascii")
wl = dat['wavelen'].values
fl = dat['flux'].values
efl = dat['flux_unc'].values
spec = Spectrum1D(
        spectral_axis=wl * u.AA, flux=fl*u.Unit('erg cm-2 s-1 AA-1'),
        uncertainty=StdDevUncertainty(efl))
spec1_gsmooth = gaussian_smooth(spec, stddev=3)

fig,ax = plt.subplots(1,1,figsize=(9,5))

z = 0.256
plt.step(wl/(1+z),fl,lw=0.3,c='k',alpha=0.5)
plt.step(spec1_gsmooth.spectral_axis/(1+z), spec1_gsmooth.flux, lw=0.5, c='k')

print(3515/(1+z))
print(6730/(1+z))
plt.axvline(x=3515/(1+z), lw=1, c='red', ls=':')
plt.axvline(x=6730/(1+z), lw=1, c='red', ls=':')

plt.ylim(0.1E-17,1.1E-17)
plt.xlim(3100/(1+z),9600/(1+z))
plt.xlabel("Rest Wavelength (AA)")
plt.ylabel("Flux Density")
plt.show()
#plt.savefig("spec_unbinned_rest.png", dpi=300)
#plt.close()
