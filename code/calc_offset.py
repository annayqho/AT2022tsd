""" Calculate the host galaxy offset """

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18
import vals

# Host coordinates from PS1
hostra = 50.045101247882215
hostdec = 8.74916033310639
hostpos = SkyCoord(hostra, hostdec, unit='deg')

# Transient coordinates from VLA
tpos = SkyCoord("03h20m10.873s", "08d44m55.739s", frame='icrs')

# Offset in arcsec
offset_arcsec = hostpos.separation(tpos).arcsec
print(offset_arcsec)

# Convert to radians
offset_rad = (offset_arcsec/3600)*np.pi/180

# Offset in physical units
dA = Planck18.angular_diameter_distance(z=vals.z).value
offset_kpc = offset_rad * dA * 1000
print(offset_kpc)
