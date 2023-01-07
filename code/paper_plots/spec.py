""" Plot the Keck/LRIS spectra """

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
import matplotlib.gridspec as gridspec
from get_spec import *
import vals
from fit_spectrum import *


def main_spec(ax, wl, flam):
    """ The main spectrum panel """

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

    # Formatting, labeling
    ax.set_yscale('log')
    plt.yticks([])
    plt.minorticks_off()
    ax.set_ylim(1.27, 3)
    ax.set_xlim(2647, 8154)
    ax.set_xlabel("Rest Wavelength ($\AA$) at $z=%s$" %vals.z)
    ax.set_ylabel("Flux (arbitrary units)")#[erg/s/cm${}^2/\AA$]")


def panels(ax, wl_range, wl, flam):
    """ Smaller panels with zoom-ins """
    choose = np.logical_and(wl>wl_range[0], wl<wl_range[1])
    x = wl[choose]
    y = flam[choose]/1E-16
    ax.step(x, y, where='mid', c='k', lw=0.5)
    plt.yticks([])
    ax.set_ylim(min(y), max(y)+0.15)
    ax.set_xlim(min(x), max(x))


def plot_lines(ax, species, col, lw=1):
    for l in wl_lines[species]:
        ax.axvline(l, lw=lw, ymin=0.9, ymax=1, color=col)


if __name__=="__main__":
    # Get data
    wl, flam, eflam = load_spec_1()
    wl = wl/(1+vals.z)
    wl_lines = get_rest_wl()

    # Initialize figure
    fig = plt.figure(figsize=(6,5))
    gs = gridspec.GridSpec(2, 3, height_ratios=(1,2))

    ax = fig.add_subplot(gs[1, :])
    main_spec(ax, wl, flam)
    ax.axvspan(3700, 3760, alpha=0.2, color='grey', lw=0)
    ax.axvspan(4800, 5050, alpha=0.2, color='grey', lw=0)
    ax.axvspan(6510, 6770, alpha=0.2, color='grey', lw=0)

    # Zoom-in of the OII (left-most) doublet
    ax = fig.add_subplot(gs[0, 0])
    panels(ax, [3700, 3760], wl, flam)
    plot_lines(ax, 'oii', vals.gc)
    ax.text(0.05, 0.95, '[O II]', ha='left', va='top', transform=ax.transAxes,
            color=vals.gc)

    # Zoom-in of the middle lines
    ax = fig.add_subplot(gs[0, 1])
    panels(ax, [4800, 5050], wl, flam)
    plot_lines(ax, 'oiii', vals.gc)
    ax.text(0.35, 0.95, '[O III]', ha='left', va='top', transform=ax.transAxes,
            color=vals.gc)
    plot_lines(ax, 'hb', vals.rc, lw=2)
    ax.text(0.00, 0.95, r'[H$\beta$]', ha='left', va='top', 
            transform=ax.transAxes, color=vals.rc)

    # Zoom-in of the right lines
    ax = fig.add_subplot(gs[0, 2])
    panels(ax, [6510, 6770], wl, flam)
    plot_lines(ax, 'ha', vals.rc, lw=2)
    ax.text(0.01, 0.9, r'[H$\alpha$]', ha='left', va='top', 
            transform=ax.transAxes, color=vals.rc)
    plot_lines(ax, 'nii', vals.gc)
    ax.text(0.3, 0.95, r'[N II]', ha='left', va='top', 
            transform=ax.transAxes, color=vals.gc)
    plot_lines(ax, 'sii', 'grey', lw=3)
    ax.text(0.8, 0.85, r'[S II]', ha='center', va='top', 
            transform=ax.transAxes, color='grey')

    #plt.tight_layout()
    #plt.show()
    plt.savefig("spec.png", dpi=200, bbox_inches='tight', pad_inches=0.1)
    plt.close()
