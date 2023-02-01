""" Calculate the minutes of exposures and number of different nights
using different facilities, after the initial IMACS detection """

from get_opt import *

# Get the number of minutes of exposures
# Can't do this yet...need Dan's updated document that includes exposure times
tel,mjd,filt,mag,emag,limmag,flare = get_flares()
jd,exp,filt,mag,emag,fujy,efujy = get_ipac()
