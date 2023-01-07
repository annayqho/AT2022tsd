""" Plot the Keck/LRIS spectra """

import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
import matplotlib.gridspec as gridspec
from get_spec import *
import vals
from fit_spectrum import *


def main_spec(ax):
    """ The main spectrum panel """
    ax = fig.add_subplot(gs[1, :])

    # Get data
    wl, flam, eflam = load_spec_1()
    wl = wl/(1+vals.z)

    # Plot
    offset = 2
    ax.step(wl, flam/1E-16+offset, where='mid', c='k', lw=0.5)#, alpha=0.5)
    ax.text(8000, offset-0.1, 
            'Keck+LRIS Sept. 23 ($\Delta t_\mathrm{rest}=13\,$d)', 
            ha='right', va='top')

    # Get second epoch
    wl, flam, eflam = load_spec_2()
    wl = wl/(1+vals.z)

    # Plot
    offset = 1.4
    ax.step(wl, flam/1E-16+offset, where='mid', c='k', lw=0.5)
    ax.text(8000, offset-0.03, 
            'Keck+LRIS Oct. 6 ($\Delta t_\mathrm{rest}=23\,$d)', 
            ha='right', va='top')


def panels(ax):
    """ Smaller panels with zoom-ins """
    # Zoom-in of the left-most line
    ax = fig.add_subplot(gs[0, 0])
    choose = np.logical_and(wl>3700, wl < 3760)
    offset = 1.4
    ax.step(wl[choose], flam[choose]/1E-16+offset, where='mid', c='k', lw=0.5)
    plt.yticks([])

    # Zoom-in of the middle lines
    ax = fig.add_subplot(gs[0, 1])
    choose = np.logical_and(wl>4800, wl < 5050)
    offset = 1.4
    ax.step(wl[choose], flam[choose]/1E-16+offset, where='mid', c='k', lw=0.5)
    plt.yticks([])

    # Zoom-in of the right lines
    ax = fig.add_subplot(gs[0, 2])
    choose = np.logical_and(wl>6530, wl < 6750)
    offset = 1.4
    ax.step(wl[choose], flam[choose]/1E-16+offset, where='mid', c='k', lw=0.5)
    plt.yticks([])


if __name__=="__main__":
    # Initialize figure
    fig = plt.figure(figsize=(6,5))
    gs = gridspec.GridSpec(2, 3, height_ratios=(1,2))


    # Mark lines
    #lines = get_rest_wl()
    #for key,line in lines.items():
    #    for l in line:
    #        ax.axvline(x=l, lw=1, alpha=0.3, c='grey')


    wl_lines['ha'] = [6562.819]
    wl_lines['hb'] = [4861.333]
    wl_lines['oii'] = [3726.032, 3728.815]
    wl_lines['oiii'] = [4958.911, 5006.843]
    wl_lines['nii'] = [6548.050, 6583.460]
    # only fit the first line...the second one has nans in it
    wl_lines['sii'] = [6716.440]#, 6730.810]

    # Formatting, labeling
    ax.set_yscale('log')
    plt.yticks([])
    plt.minorticks_off()
    ax.set_ylim(1.27, 3)
    ax.set_xlim(2647, 8154)
    ax.set_xlabel("Rest Wavelength ($\AA$) at $z=%s$" %vals.z)
    ax.set_ylabel("Flux (arbitrary units)")#[erg/s/cm${}^2/\AA$]")

    plt.tight_layout()
    plt.show()
