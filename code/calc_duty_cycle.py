"""
Calculate the duty cycle of the flares for different limiting magnitudes
"""

from get_opt import *

# Get the optical photometry
dat = get_full_opt()

# I think it's reasonably convincing that the flaring starts with the first
# ZTF detection (aka the first flare detection).

# I think we can say that the flaring "ends" with the last flare detection.

# Get the 

# For a given limiting magnitude, sum the number of exposures

# Get the number of minutes of exposures
# Can't do this yet...need Dan's updated document that includes exposure times
tel,mjd,filt,mag,emag,limmag,flare = get_flares()
jd,exp,filt,mag,emag,fujy,efujy = get_ipac()
