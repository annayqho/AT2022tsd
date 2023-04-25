""" Figure showing the position of the transient in its host galaxy,
as well as the position of the host galaxy in the M*-SFR plane """

import numpy as np
import matplotlib.pyplot as plt
import astropy.wcs
from astropy import coordinates as coords
from astropy.io import fits
from astropy.visualization import make_lupton_rgb
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
import vals

ddir = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/host/"
tcol = vals.lgrb_col
tcol = 'k'

def get_host_phot(imsize):
    # Position of the transient in the host galaxy
    ra = vals.ra
    dec = vals.dec

    # Figure out pos from header
    rim = fits.open(ddir + "cutout_rings.v3.skycell.1513.011.stk.r.unconv.fits")
    head = rim[0].header
    wcs = astropy.wcs.WCS(head)
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

    return gcut,rcut,icut,zcut


if __name__=="__main__":
    # Initiate figure
    fig,axarr = plt.subplots(1,2,figsize=(8,4))

    ax = axarr[0]

    # Plot host galaxy
    imsize = 50
    g,r,i,z = get_host_phot(imsize)
    ax.imshow(r, origin='lower', cmap='Greys', vmin=-200, vmax=500)

    # Mark position of transient
    ax.plot([imsize, imsize], [imsize, imsize-8], c=tcol, lw=1)
    ax.plot([imsize, imsize+8], [imsize, imsize], c=tcol, lw=1)
    ax.text(imsize+1, imsize-2, 'AT2022tsd', ha='left', va='top', fontsize=10)

    # Mark compass
    imsize = 100
    ax.plot((imsize-10,imsize-10), (imsize-10,imsize-20), color='k', lw=2)
    ax.text(
            imsize-10, imsize-23, "S", color='k', fontsize=16,
            horizontalalignment='center', verticalalignment='top')
    ax.plot((imsize-10,imsize-20), (imsize-10,imsize-10), color='k', lw=2)
    ax.text(
            imsize-23, imsize-10, "E", color='k', fontsize=16,
            horizontalalignment='right', verticalalignment='center')
    ax.axis('off')

    # Mark image scale : 0.25 arcsec per pixel
    x = 7
    y = 10
    x2 = x + 5/0.25
    ax.plot((x,x2), (y,y), c='k', lw=2)
    ax.text((x2+x)/2, y*1.1, "5''", color='k', fontsize=16,
            verticalalignment='bottom', horizontalalignment='center')
    ax.text((x2+x)/2, y/1.1, "(32 kpc)", color='k', fontsize=16,
            verticalalignment='top', horizontalalignment='center')

    #ax = axarr[1]

    plt.show()
