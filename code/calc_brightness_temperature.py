""" Calculate the brightness temperature of the ULTRASPEC flares """
import vals
from get_opt import *

# Constants
c = 3E10
k = 1.38E-16
Gamma = 1

# Source size
dt = 30 # observer frame
dR = c * dt * Gamma**2
dTheta = dR / vals.dL_cm
print(dTheta * (180 / np.pi) * 3600 * 1E6)

# Intensity
dat = get_full_opt()
tel = dat['#instrument'].values
flux = dat['flux_extcorr'].values
isflare = dat['isflare'].values
filt = dat['flt'].values

choose = np.logical_and(tel=='TNT/ULTRASPEC', isflare)
peak_ind = np.argmax(flux[choose])
peak_filt = filt[choose][peak_ind]
Snu = flux[choose][peak_ind]*1E-6*(1+vals.z) # Jy
Inu = Snu * 1E-23 / (np.pi * dTheta**2)

# Frequency
nu = 3E18/(vals.sdss_pivot[peak_filt])

# Brightness temperature
TB = (Inu*c**2) / (2*k*nu**2)
print(TB/1E10)
