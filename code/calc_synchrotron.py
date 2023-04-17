""" Synchrotron calculations 

Following Fender et al. (1997) for GRS 1915+105
and Readhead synchrotron notes
"""

import numpy as np

# Input: Fender values
# R = 9E12
# L = 1E36
# alpha = 0
# nu1 = 1E14
# nu2 = 15*1E9
# obs_nu = 1.4E14

# Input: our values
R = 1E12
L = 1E44
alpha = -1.6
nu1 = 1E15
nu2 = 1E12
obs_nu = 1E15

# Constants
A = 1.586E12

# Radius: 1 light-minute, which is 10^12 cm
V = (4/3) * np.pi * R**3

# Values of g...uncertain here
g = ((2*alpha+2)/(2*alpha+1)) * ((nu2**(alpha+1/2)-nu1**(alpha+1/2))/(nu2**(alpha+1)-nu1**(alpha+1)))

# Equiparititon magnetic field
B = (8*np.pi*A*g*L/V)**(2/7)
print(B)

# Eqiuparittion energy
U = V*B**2/(8*np.pi) + A*g*L*B**(-3/2)
print(U)

# For a given magnetic field strength, what is the required particle energy
# at optical frequencies, say 1E15 Hz?
nu_g = (2.8 / B)*1E6 # nu_g = 2.8 MHz / Gauss; this value is in Hz
gamma_e = np.sqrt((obs_nu/nu_g))
particle_energy_erg = gamma_e * 9E-28 * (3E10)**2
particle_energy_eV = particle_energy_erg*6.242E11
print(particle_energy_eV/1E9)

tcool = 6*np.pi*9E-28*3E10 / (6.65E-25 * B**2 * gamma_e)
print(tcool) # cooling time in hours
