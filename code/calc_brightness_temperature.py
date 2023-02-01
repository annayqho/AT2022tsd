""" Calculate the brightness temperature of the ULTRASPEC flares """
import vals
from get_opt import *

# Constants
c = 3E10
k = 1.38E-16

# Source size
dt = 30/(1+vals.z) # 30 seconds in the rest frame
dR = c * dt
dTheta = dR / vals.dL_cm
print(dTheta * (180 / np.pi) * 3600 * 1E6)

# Intensity
tel,mjd,filt,mag,emag,limmag,flare = get_flares()
choose = np.logical_and.reduce((emag<99, tel=='ULTRASPEC', ~np.isnan(mag)))
peak_ind = np.argmin(mag[choose])
peak_filt = filt[choose][peak_ind]
peak_mag = mag[choose][peak_ind]
ext_corr_peak = peak_mag-vals.ext[peak_filt]
Snu = 10**((ext_corr_peak-8.90)/(-2.5)) # Jy
Inu = Snu * 1E-23 / (np.pi * dTheta**2)

# Frequency
nu = 3E18/(vals.sdss_pivot[peak_filt])

# Brightness temperature
TB = (Inu*c**2) / (2*k*nu**2)
print(TB/1E10)
