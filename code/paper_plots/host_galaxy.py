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

    return [gcut,rcut,icut]#,zcut


def get_host_phot_lris(imsize):
    """ Get host photometry from LRIS """

    # Basic parameters
    ddir = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/opt/LRIS/imaging/"
    ra = vals.ra
    dec = vals.dec

    # Image data that we will use
    ifits = "lrisr20230117_ZTF22abftjko_i.fits" # i-band image
    gfits = "lrisb20230117_ZTF22abftjko_G.fits" # g-band image
    iim = fits.open(ddir+ifits) # i-band image
    gim = fits.open(ddir+gfits) # g-band image
    ihead = iim[0].header # i-band header
    ghead = gim[0].header # g-band header
    iproj = WCS(ihead) # i-band projection
    gproj = WCS(ghead) # g-band projection
    idat = iim[0].data # i-band data
    gdat = gim[0].data # g-band data

    # Extract the i-band image
    target_pix = iproj.all_world2pix([(np.array([ra,dec], np.float_))], 1)[0]
    xpos = target_pix[0]
    ypos = target_pix[1]
    icut = idat[int(ypos-imsize):int(ypos+imsize),
                int(xpos-imsize):int(xpos+imsize)]

    # Resample the g-band image onto the i-band image
    greproj, _ = reproject_adaptive((gdat, gproj),
                        iproj, shape_out=idat.shape)
    gcut = greproj[int(ypos-imsize):int(ypos+imsize),
                int(xpos-imsize):int(xpos+imsize)]

    return icut,gcut,gcut


if __name__=="__main__":
    # Initiate figure
    fig,axarr = plt.subplots(1,2,figsize=(8,4))

    ax = axarr[0]

    # Plot host galaxy
    imsize = 50
    r,g,b = get_host_phot_lris(imsize) # r:0-1000; g:0-500
    #r,g,b = get_host_phot_ps1(imsize) # r:0-1000; g:0-500
    rgb = make_lupton_rgb(r/2, g, b, Q=2, stretch=0.1)
    ax.imshow(rgb, origin='lower', vmin=100, vmax=500)#, cmap='Greys', vmin=-200, vmax=500)

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

    # Plot the Taggart sample (scraped from 20xnd paper)
    ax = axarr[1]
    dat = np.loadtxt("taggart_sample.txt", delimiter=',')
    ax.scatter(dat[:,0], dat[:,1], marker='>', c=vals.sn_col, s=10, label='CCSN')
    dat = np.loadtxt("taggart_lgrb.txt", delimiter=',')
    ax.scatter(dat[:,0], dat[:,1], marker='+', c=vals.lgrb_col, s=40, label='LGRB')

    # Plot LFBOTs
    x = [1E8, 1E8, 1.3E7, 1.4E9, 3.1E8, 10**9.31]
    y = [6.93E-3, 2E-2, 4E-3, 2.2E-1, 6.2, 3.34]
    ax.scatter(x, y, marker='D', facecolor=vals.cow_col, edgecolor='k', label='LFBOT') 
    ax.text(x[-1]*1.1, y[-1], 'AT2022tsd', c=vals.cow_col, ha='left', 
            va='bottom', fontsize=8, fontweight='bold')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r"Stellar Mass ($M_\odot$)")
    ax.set_ylabel(r"Star formation rate ($M_\odot\,\mathrm{yr}^{-1}$)")
    ax.legend(loc='lower right', fontsize=8)
    plt.tight_layout()
    plt.show()
    #plt.savefig("host_galaxy.png", dpi=300)
    #plt.close()
