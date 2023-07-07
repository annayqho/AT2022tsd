""" Fit a blackbody to the peak-light SED """

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import vals


def bb_func(wl_AA,T,R):
    """
    blackbody function

    Parameters
    ----------
    wl: wavelength in angstroms
    T: temperature in Kelvin
    R: radius in cm

    Returns
    -------
    flux density in units of uJy
    """
    d = vals.dA_cm
    h = 6.626E-27
    c = 3E10
    k = 1.38E-16
    
    # Convert wl to cm
    wl = wl_AA / 1E8

    Blam = (2*h*c**2/wl**5) * (1/(np.exp(h*c/(wl*k*T))-1))
    flam = Blam * np.pi * R**2 / d**2
    # in units of uJy
    fnu = wl**2 * flam / c
    return fnu / 1E-23 / 1E-6


# g, r, i: wl (AA), flux (uJy), eflux (uJy)
wl_obs = np.array([vals.ztf_leff['g'], vals.ztf_leff['r'], vals.ps1_leff['i']])
# Convert to rest
wl = wl_obs / (1+vals.z)
f = np.array([70, 46, 36])
ef = np.array([6, 5, 6])

plt.errorbar(wl,f,ef,fmt='o')

# Fit the blackbody
xmod = np.linspace(3000,7000)
ymod = bb_func(xmod, 40000, 1E14)
popt, pcov = curve_fit(
        bb_func, wl, f, sigma=ef, absolute_sigma=True, p0=[30000,1E14])
plt.plot(xmod, bb_func(xmod, *popt))
plt.show()
