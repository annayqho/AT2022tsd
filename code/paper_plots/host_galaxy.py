""" Figure showing the position of the transient in its host galaxy,
as well as the position of the host galaxy in the M*-SFR plane """

import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from reproject import reproject_interp, reproject_adaptive
from astropy import coordinates as coords
from astropy.io import fits
from astropy.visualization import make_lupton_rgb
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
import vals

tcol = vals.lgrb_col
tcol = 'k'

def get_host_phot_ps1(imsize):
    """ Get host photometry from PS1 """
    ddir = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/host/ps1/"
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
    ddir = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/opt/LRIS/imaging/"
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


if __name__=="__main__":
    # Initiate figure
    fig,axarr = plt.subplots(1,2,figsize=(8,4))

    ax = axarr[0]

    # Plot host galaxy
    imsize = 50
    i,g,u = get_host_phot_lris(imsize)
    r = (i+g)/2  # max: 1000
    # max for g: 600, 650
    b = (g*1.5+u)/2  # max for b: 500, 550

    #r,g,b = get_host_phot_ps1(imsize) # r:0-1000; g:0-500
    #rgb = make_lupton_rgb(r/1.8, g/1.1, b, stretch=100, Q=5, minimum=13)

    # Try gri
    #rgb = make_lupton_rgb(i/2.5, r/1.9, b, stretch=100, Q=4, minimum=10)
    rgb = make_lupton_rgb(i/2.5, r/1.9, b, stretch=100, Q=5, minimum=10)
    ax.imshow(rgb, origin='lower')

    markcol = 'white'

    # Mark position of transient
    ax.plot([imsize, imsize], [imsize, imsize-8], c=markcol, lw=1)
    ax.plot([imsize, imsize+8], [imsize, imsize], c=markcol, lw=1)
    ax.text(imsize+1, imsize-2, 'AT2022tsd', ha='left', va='top', fontsize=10,
            color=markcol)

    # Mark compass
    imsize = 100
    ax.plot((imsize-10,imsize-10), (imsize-10,imsize-20), color=markcol, lw=2)
    ax.text(
            imsize-10, imsize-23, "S", color=markcol, fontsize=16,
            horizontalalignment='center', verticalalignment='top')
    ax.plot((imsize-10,imsize-20), (imsize-10,imsize-10), color=markcol, lw=2)
    ax.text(
            imsize-23, imsize-10, "E", color=markcol, fontsize=16,
            horizontalalignment='right', verticalalignment='center')
    ax.axis('off')

    ax.text(0.01, 0.99,"Keck/LRIS $ugI$",fontsize=15,transform=ax.transAxes,
        horizontalalignment='left', va='top', color='white')

    # Mark image scale : 0.25 arcsec per pixel
    x = 7
    y = 10
    x2 = x + 5/0.25
    ax.plot((x,x2), (y,y), color=markcol, lw=2)
    ax.text((x2+x)/2, y*1.1, "5''", color=markcol, fontsize=16,
            verticalalignment='bottom', horizontalalignment='center')
    ax.text((x2+x)/2, y/1.1, "(21 kpc)", color=markcol, fontsize=16,
            verticalalignment='top', horizontalalignment='center')

    # Plot the Taggart sample (scraped from 20xnd paper)
    ax = axarr[1]
    dat = np.loadtxt("taggart_sample.txt", delimiter=',')
    ax.scatter(dat[:,0], dat[:,1], marker='>', c=vals.sn_col, s=10, label='CCSN')
    dat = np.loadtxt("taggart_lgrb.txt", delimiter=',')
    ax.scatter(dat[:,0], dat[:,1], marker='+', c=vals.lgrb_col, s=40, label='LGRB')

    # Add 18cow
    x = 1.42E9
    xbot = (1.42-0.29)*1E9
    xtop = (1.42+0.17)*1E9
    y = 0.22
    ybot = 0.22-0.04
    ytop = 0.22+0.03
    ax.scatter(x, y, marker='D', facecolor=vals.cow_col, edgecolor='k', label='LFBOT', zorder=100)
    ax.hlines(y, xbot, xtop, color=vals.cow_col)
    ax.vlines(x, ybot, ytop, color=vals.cow_col)

    # Add CSS161010
    x = 2E7
    xbot = 1E7
    xtop = 3E7
    y = 1E-2
    ybot = 0.3E-2
    ytop = 2E-2
    ax.scatter(x, y, marker='D', facecolor=vals.cow_col, edgecolor='k', label=None, zorder=100)
    ax.hlines(y, xbot, xtop, color=vals.cow_col)
    ax.vlines(x, ybot, ytop, color=vals.cow_col)

    # Add Koala
    x = 5.1E8
    xbot = (5.1-2.0)*1E8
    xtop = (5.1+3.4)*1E8
    y = 6.8
    ybot = (6.8-4.6)
    ytop = (6.8+3.7)
    ax.scatter(x, y, marker='D', facecolor=vals.cow_col, edgecolor='k', label=None, zorder=100)
    ax.hlines(y, xbot, xtop, color=vals.cow_col)
    ax.vlines(x, ybot, ytop, color=vals.cow_col)

    # Add AT2020xnd
    x = 8E7
    xbot = 3E7
    xtop = 3E8
    y = 0.02
    ybot = 0.02-0.005
    ytop = 0.02+0.005
    ax.scatter(x, y, marker='D', facecolor=vals.cow_col, edgecolor='k', label=None, zorder=100)
    ax.hlines(y, xbot, xtop, color=vals.cow_col)
    ax.vlines(x, ybot, ytop, color=vals.cow_col)

    # Add 2020mrf
    x = 10**7.94
    xbot = 10**(7.94-0.39)
    xtop = 10**(7.94+0.22)
    y = 6.93*1E-3
    ybot = (6.93-0.27)*1E-3
    ytop = (6.93+3.90)*1E-3
    ax.scatter(x, y, marker='D', facecolor=vals.cow_col, edgecolor='k', label=None, zorder=100)
    ax.hlines(y, xbot, xtop, color=vals.cow_col)
    ax.vlines(x, ybot, ytop, color=vals.cow_col)

    # Add AT2022tsd
    x = 10**(9.96)
    xtop = 10**(9.96+0.06)
    xbot = 10**(9.96-0.09)
    y = 0.55
    ytop = 0.55+1.36
    ybot = 0.55-0.19
    ax.scatter(x, y, marker='D', facecolor=vals.cow_col, edgecolor='k', label=None, zorder=100)
    ax.hlines(y, xbot, xtop, color=vals.cow_col)
    ax.vlines(x, ybot, ytop, color=vals.cow_col)
    ax.text(x*1.1, y/1.1, 'AT2022tsd', c=vals.cow_col, ha='left',
            va='top', fontsize=8, fontweight='bold')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r"Stellar Mass ($M_\odot$)")
    ax.set_ylabel(r"Star formation rate ($M_\odot\,\mathrm{yr}^{-1}$)")
    ax.legend(loc='lower right', fontsize=8)
    plt.tight_layout()
    #plt.show()
    plt.savefig("host_galaxy.png", dpi=300)
    plt.close()
