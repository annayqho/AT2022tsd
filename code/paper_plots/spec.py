""" Plot the Keck/LRIS spectra """

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.optimize import curve_fit
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
import matplotlib.gridspec as gridspec
from get_spec import *
import vals
from fit_spectrum import *

single_width = 40
multi_width = 130


def gauss(x, sigma, A, b):
    """ Gaussian distribution """
    return A*np.exp(-(x-b)**2/(2*sigma**2))


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
    ax.set_xlabel("$\lambda_\mathrm{rest}$ ($\AA$) at $z=%s$" %vals.z)
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
    wl_lines = get_rest_wl()
    for l in wl_lines[species]:
        ax.axvline(l, lw=lw, ymin=0.9, ymax=1, color=col)


def fig_for_paper():
    # Get data
    wl, flam, eflam = load_spec_1()
    wl = wl/(1+vals.z)
    wl2, flam2, eflam2 = load_spec_2()
    wl2 = wl2/(1+vals.z)
    wl_lines = get_rest_wl()

    # Initialize figure
    fig = plt.figure(figsize=(6,7))
    gs = gridspec.GridSpec(3, 3, height_ratios=(1,2,1))

    ax = fig.add_subplot(gs[1, :])
    main_spec(ax, wl, flam)
    ax.axvspan(3700, 3760, alpha=0.2, color='grey', lw=0)
    ax.axvspan(4800, 5050, alpha=0.2, color='grey', lw=0)
    ax.axvspan(6510, 6770, alpha=0.2, color='grey', lw=0)

    # Zoom-in of the OII (left-most) doublet
    ax = fig.add_subplot(gs[0, 0])
    panels(ax, [3725-single_width, 3725+single_width], wl, flam+2E-17)
    panels(ax, [3725-single_width, 3725+single_width], wl2, flam2)
    plot_lines(ax, 'oii', vals.gc)
    ax.text(0.05, 0.95, '[O II]', ha='left', va='top', transform=ax.transAxes,
            color=vals.gc)
    ax.set_ylim(0, 1)

    # Zoom-in of the middle lines
    ax = fig.add_subplot(gs[0, 1])
    panels(ax, [4930-multi_width, 4930+multi_width], wl, flam+2E-17)
    panels(ax, [4930-multi_width, 4930+multi_width], wl2, flam2)
    plot_lines(ax, 'oiii', vals.gc)
    ax.text(0.35, 0.95, '[O III]', ha='left', va='top', 
            transform=ax.transAxes, color=vals.gc)
    plot_lines(ax, 'hb', vals.rc, lw=2)
    ax.text(0.00, 0.95, r'[H$\beta$]', ha='left', va='top', 
            transform=ax.transAxes, color=vals.rc)
    ax.set_ylim(0, 0.8)

    # Zoom-in of the right lines
    ax = fig.add_subplot(gs[0, 2])
    panels(ax, [6640-multi_width, 6640+multi_width], wl, flam+2E-17)
    panels(ax, [6640-multi_width, 6640+multi_width], wl2, flam2)
    plot_lines(ax, 'ha', vals.rc, lw=2)
    ax.text(0.01, 0.9, r'[H$\alpha$]', ha='left', va='top', 
            transform=ax.transAxes, color=vals.rc)
    plot_lines(ax, 'nii', vals.gc)
    ax.text(0.3, 0.95, r'[N II]', ha='left', va='top', 
            transform=ax.transAxes, color=vals.gc)
    plot_lines(ax, 'sii', 'grey', lw=3)
    ax.text(6725, 0.85, r'[S II]', ha='center', va='top', color='grey')
    ax.set_ylim(0, 1.2)

    # Get data assuming z=0
    wl, flam, eflam = load_spec_1()
    wl2, flam2, eflam2 = load_spec_2()
    wl_lines = get_rest_wl()

    # Zoom-in of H-alpha at z=0 (6563 AA)
    ax = fig.add_subplot(gs[2, 0])
    panels(ax, [6560-single_width, 6560+single_width], wl, flam+0.7E-17)
    panels(ax, [6560-single_width, 6560+single_width], wl2, flam2)
    plot_lines(ax, 'ha', vals.rc)
    ax.text(0.05, 0.95, r'H$\alpha$', ha='left', va='top', 
            transform=ax.transAxes, color=vals.rc)
    #ax.set_xticks([6520, 6550, 6580])
    #ax.set_xticklabels([6520, 6550, 6580])
    ax.set_ylim(0, 0.35)
    ax.set_xlabel("$\lambda_\mathrm{obs}$ ($\AA$)")

    # Zoom-in of He I 5875
    ax = fig.add_subplot(gs[2, 1])
    panels(ax, [5875-single_width, 5875+single_width], wl, flam+2E-17)
    panels(ax, [5875-single_width, 5875+single_width], wl2, flam2)
    plot_lines(ax, 'hei', vals.rc)
    ax.text(0.05, 0.95, r'HeI 5875', ha='left', va='top', 
            transform=ax.transAxes, color=vals.rc)
    ax.set_ylim(0, 0.6)
    ax.set_xlabel("$\lambda_\mathrm{obs}$ ($\AA$)")
    #ax.set_xticks([5400, 5700, 6200])
    #ax.set_xticklabels([5400, 5700, 6200])

    #plt.tight_layout()

    # Zoom-in of He II 4686
    ax = fig.add_subplot(gs[2, 2])
    panels(ax, [4686-single_width, 4686+single_width], wl, flam+2E-17)
    panels(ax, [4686-single_width, 4686+single_width], wl2, flam2)
    plot_lines(ax, 'heii', vals.rc)
    ax.text(0.05, 0.95, r'HeII 4686', ha='left', va='top', 
            transform=ax.transAxes, color=vals.rc)
    ax.set_ylim(0, 1)
    ax.set_xlabel("$\lambda_\mathrm{obs}$ ($\AA$)")
    #ax.set_xticks([4660, 4680, 4700, 4720])
    #ax.set_xticklabels([4660, 4680, 4700, 4720])

    plt.tight_layout()
    plt.show()
    #plt.savefig("spec.png", dpi=200, bbox_inches='tight', pad_inches=0.1)
    #plt.close()


if __name__=="__main__":
    fig_for_paper()
    # measure the offset in position from the HeII feature
    # fig,ax = plt.subplots(1,1)
    # #wl, flam, eflam = load_spec_1()
    # wl, flam, eflam = load_spec_2()
    # # normalize and plot
    # choose = np.logical_and(wl>4690, wl < 4720)
    # med = np.median(flam[choose])
    # flam = flam/med
    # choose = np.logical_and(wl>4660, wl < 4700)
    # ax.step(wl[choose], flam[choose]-1, where='mid', lw=0.5, c='grey')
    # # fit a Gaussian
    # popt, pcov = curve_fit(gauss, wl[choose], flam[choose]-1, p0=(5, 4, 4700))
    # xvals = np.linspace(4665, 4700)
    # yvals = gauss(xvals, *popt)
    # ax.plot(xvals, yvals, c='k', lw=2)
    # # peak of spec1: 4683.54
    # # peak of spec2: 4686.68

    # w = 4686
    # 
    


