""" Calculate the extinction-corrected flare color based on 
Keck/LRIS data """

import numpy as np
import vals

def get_mags():
    """ get extinction corrected magnitudes """
    u = 25.10-vals.ext['u']
    eu = 0.30
    i = 22.99-vals.ext['i']
    ei = 0.07
    return u, i, eu, ei

def calc_color_mag():
    u, i, eu, ei = get_mags()
    ui = u-i
    eui = np.sqrt(eu**2+ei**2)
    return ui, eui

def calc_color_spindex():
    """ Calculate the color as the spectral index, fnu ~ nu^{beta} """
    nu_u = 3E18 / vals.keck_leff['u'] 
    nu_I = 3E18 / vals.keck_leff['I'] 
    ui, eui = calc_color_mag()
    # Convert to a flux ratio
    fratio = 10**(ui/(-2.5))
    # Solve for beta
    nuratio = nu_u / nu_I
    beta = np.log(fratio)/np.log(nuratio)
    # Equivalent expression for beta
    beta = (1/np.log10(nuratio)) * (ui)/(-2.5)
    # I think the uncertainty is: sqrt((em1/2.5)**2 + (em2/2.5)**2)
    u, i, eu, ei = get_mags()
    ebeta = np.sqrt((ei/2.5)**2+(eu/2.5)**2)
    print(beta, ebeta)

if __name__=="__main__":
    calc_color_spindex()
    


