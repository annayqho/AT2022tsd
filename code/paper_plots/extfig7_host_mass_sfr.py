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
sys.path.append("..")
import vals
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 7 # The maximum allowed for ED figures

tcol = vals.lgrb_col
tcol = 'k'


if __name__=="__main__":
    figwidth_mm = 89 # Nature standard
    figwidth_in = (figwidth_mm/10)/2.54 # in inches

    # Initiate figure
    fig,ax = plt.subplots(1,1,figsize=(figwidth_in,figwidth_in))

    # Plot the Taggart sample (scraped from 20xnd paper)
    dat = np.loadtxt("../../data/taggart_sample.txt", delimiter=',')
    ax.scatter(dat[:,0], dat[:,1], marker='>', c=vals.sn_col, s=10, label='CCSN')
    dat = np.loadtxt("../../data/taggart_lgrb.txt", delimiter=',')
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
            va='top', fontweight='bold')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r"Stellar Mass ($M_\odot$)")
    ax.set_ylabel(r"Star formation rate ($M_\odot$ yr${}^{-1}$)")
    ax.legend(loc='lower right')
    plt.tight_layout()
    #plt.show()
    plt.savefig("host_galaxy.eps", dpi=300, bbox_inches='tight', pad_inches=0)
    plt.close()
