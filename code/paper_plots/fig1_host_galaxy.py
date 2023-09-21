""" Figure showing the position of the transient in its host galaxy,
as well as the position of the host galaxy in the M*-SFR plane """

import numpy as np
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['pdf.fonttype']=42
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from reproject import reproject_interp, reproject_adaptive
from astropy import coordinates as coords
from astropy.io import fits
from astropy.visualization import make_lupton_rgb
import sys
sys.path.append("..")
import vals

tcol = vals.lgrb_col
tcol = 'k'

def get_host_phot_ps1(imsize):
    """ Get host photometry from PS1 """
    ddir = "../../data/host/ps1/"
    # Position of the transient in the host galaxy
    ra = vals.ra
    dec = vals.dec

    # Figure out pos from header
    rim = fits.open(ddir + "cutout_rings.v3.skycell.1513.011.stk.r.unconv.fits")
    head = rim[0].header
    wcs = WCS(head)
    target_pix = wcs.all_world2pix([(np.array([ra,dec], np.float_))], 1)[0]
    xpos = target_pix[0]
    ypos = target_pix[1]

    # Plot transient in its host galaxy
    g = fits.open( # range 20~200 perhaps
            ddir + "cutout_rings.v3.skycell.1513.011.stk.g.unconv.fits")[0].data
    r = fits.open( # range 20-400 perhaps
            ddir + "cutout_rings.v3.skycell.1513.011.stk.r.unconv.fits")[0].data
    i = fits.open( # range 20-800 perhaps
            ddir + "cutout_rings.v3.skycell.1513.011.stk.i.unconv.fits")[0].data
    z = fits.open( # maybe similar to i
            ddir + "cutout_rings.v3.skycell.1513.011.stk.z.unconv.fits")[0].data

    gcut = g[int(ypos-imsize):int(ypos+imsize),int(xpos-imsize):int(xpos+imsize)]
    rcut = r[int(ypos-imsize):int(ypos+imsize),int(xpos-imsize):int(xpos+imsize)]
    icut = i[int(ypos-imsize):int(ypos+imsize),int(xpos-imsize):int(xpos+imsize)]
    zcut = z[int(ypos-imsize):int(ypos+imsize),int(xpos-imsize):int(xpos+imsize)]

    return [gcut,rcut,icut]#,zcut


def get_host_phot_lris(imsize):
    """ Get host photometry from LRIS """

    # Basic parameters
    ddir = "../../data/opt/LRIS/imaging/"
    ra = vals.ra
    dec = vals.dec

    # Image data that we will use
    ifits = "i.resamp.fits" # i-band image
    gfits = "g.resamp.fits" # g-band image
    ufits = "u.resamp.fits" # g-band image

    out = []
    for imf in [ifits,gfits,ufits]:
        im = fits.open(ddir+imf) 
        head = im[0].header 
        wcs = WCS(head)
        dat = im[0].data 

        # Extract the cutout
        target_pix = wcs.all_world2pix([(np.array([ra,dec], np.float_))], 1)[0]
        xpos = target_pix[0]
        ypos = target_pix[1]
        cut = dat[int(ypos-imsize):int(ypos+imsize),
                   int(xpos-imsize):int(xpos+imsize)]
        out.append(cut)

    return out[0], out[1], out[2]
