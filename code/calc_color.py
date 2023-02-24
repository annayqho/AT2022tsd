""" Calculate the extinction-corrected flare color based on 
Keck/LRIS data """

import numpy as np
from get_opt import *
import vals


def calc_color_spindex(dat, filtblue, filtred):
    """ Calculate the color as the spectral index, fnu ~ nu^{beta} """
    nu_blue = 3E18 / vals.keck_leff[filtblue] 
    nu_red = 3E18 / vals.keck_leff[filtred] 

    mag_blue = dat[dat['flt']==filtblue]['mag']-vals.ext[filtblue]
    mag_red = dat[dat['flt']==filtred]['mag']-vals.ext[filtred]
    emag_blue = dat[dat['flt']==filtblue]['emag']
    emag_red = dat[dat['flt']==filtred]['emag']

    bluered = mag_blue.values-mag_red.values
    ebluered = np.sqrt(emag_blue.values**2+emag_red.values**2)

    keep = np.logical_and(emag_blue.values<99, emag_red.values<99)
    bluered = bluered[keep]
    ebluered = ebluered[keep]

    # Convert to a flux ratio
    fratio = 10**(bluered/(-2.5))
    # Solve for beta
    nuratio = nu_blue / nu_red
    beta = np.log(fratio)/np.log(nuratio)
    # Equivalent expression for beta
    beta = (1/np.log10(nuratio)) * (bluered)/(-2.5)
    # I think the uncertainty is: sqrt((em1/2.5)**2 + (em2/2.5)**2)
    ebeta = np.sqrt((emag_red.values[keep]/2.5)**2+(emag_blue.values[keep]/2.5)**2)

    # Get times of observations
    t = dat[dat['flt']==filtblue]['mjdstart'].values-\
            dat[dat['flt']==filtblue]['mjdstart'].values[0]

    return beta, ebeta


if __name__=="__main__":
    dat = get_full_opt()
    choose = dat['#instrument']=='KeckI/LRIS'
    dat = dat[choose]

    # Colors for the flares
    choose = np.logical_and(dat['mjdstart']>59942, dat['mjdstart']<59943)
    calc_color_spindex(dat[choose], 'u', 'i')

    # Colors for the gi flares
    choose = np.logical_and(dat['mjdstart']>59871, dat['mjdstart']<59872)
    calc_color_spindex(dat[choose], 'g', 'i')

