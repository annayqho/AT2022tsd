""" Plot the Keck/LRIS spectra """

import pandas as pd
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 7 # The maximum allowed for ED figures
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.optimize import curve_fit
import sys
sys.path.append("..")
import matplotlib.gridspec as gridspec
from get_spec import *
import vals
from fit_spectrum import *

single_width = 40
multi_width = 130


def gauss(x, sigma, A, b):
    """ Gaussian distribution """
    return A*np.exp(-(x-b)**2/(2*sigma**2))


def main_spec(ax, wl, flam, eflam, wl2, flam2, eflam2):
    """ The main spectrum panel """
    labels = ['Keck+LRIS Sept. 23 ($\Delta t_\mathrm{rest}=13$ d)',
              'Keck+LRIS Oct. 6 ($\Delta t_\mathrm{rest}=23$ d)']
    label_offsets = [0.00, 0.03]
    offsets = [2, 1]
    x = [wl, wl2]
    y = [flam, flam2]
    ey = [eflam, eflam2]

    for i in np.arange(2):
        print(i)
        # Plot
        fac = max(y[i][x[i]>6000]) # scale by peak of Halpha
        yscaled = y[i]/fac+offsets[i]
        # Plot in the observer frame
        #ax.step(x[i]*(1+vals.z),yscaled,where='mid',c='lightgrey',lw=0.5)
        xb,yb,eyb = load_binned_spec(x[i]*(1+vals.z), y[i], ey[i])
        ybscaled = yb/fac+offsets[i]
        ax.step(xb, ybscaled, where='mid', c='k', lw=0.5)
        ax.text(8000*(1+vals.z), offsets[i]-label_offsets[i], labels[i],
                ha='right', va='top')

    # Formatting, labeling
    plt.yticks([])
    plt.minorticks_off()
    ax.set_ylim(0.77, 3)
    ax.set_xlim(2450*(1+vals.z), 8160*(1+vals.z))
    ax.set_xlabel("$\lambda_\mathrm{obs}$ ($\AA$)")
    ax.set_ylabel("Flux (arbitrary units)")#[erg/s/cm${}^2/\AA$]")

    # Create a secondary axis
    ax2 = ax.twiny()
    ax2.set_xlabel("$\lambda_\mathrm{rest}$ ($\AA$) at $z=%s$" %vals.z)
    x_f = lambda x_i: x_i/(1+vals.z)
    xmin, xmax = ax.get_xlim()
    ax2.set_xlim((x_f(xmin), x_f(xmax)))
    ax2.plot([],[])


def panels(ax, wl_range, wl, flam, shift):
    """ Smaller panels with zoom-ins """
    choose = np.logical_and(wl>wl_range[0], wl<wl_range[1])
    x = wl[choose]
    y = flam[choose]/max(flam[choose])
    ax.step(x, y-shift, where='mid', c='k', lw=0.5)
    plt.yticks([])
    ax.set_ylim(-1, 1.3)
    xmin, xmax = ax.get_xlim()
    ax.set_xlim(xmin+5, xmax-5)


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
    figwidth_mm = 183 # Nature standard
    figwidth_in = (figwidth_mm/10)/2.54 # in inches
    
    fig = plt.figure(figsize=(figwidth_in,figwidth_in*(7/6)))
    gs = gridspec.GridSpec(3, 3, height_ratios=(1,2,1))

    ax = fig.add_subplot(gs[1, :])
    main_spec(ax, wl, flam, eflam, wl2, flam2, eflam2)

    # Get the corners of the main box
    xlim_ax = ax.get_xlim()
    ylim_ax = ax.get_ylim()

    # Regions of the obs-frame lines
    ls = [4868, 5876, 6560]
    for l in ls:
        ax.axvspan(l-single_width, l+single_width, color='lavender', lw=0)

    # Zoom-in of the OII (left-most) doublet
    ax2 = fig.add_subplot(gs[0, 0])
    panels(ax2, [3725-single_width, 3725+single_width], wl, flam, 0)
    panels(ax2, [3725-single_width, 3725+single_width], wl2, flam2, 1)
    plot_lines(ax2, 'oii', vals.gc)
    ax2.text(0.38, 0.98, '[O II]', ha='left', va='top', 
            transform=ax2.transAxes, color=vals.gc)
    ax2.set_xlabel("$\lambda_\mathrm{rest}$ $(\AA)$")

    # Get the corners of the panel
    xlim_ax2 = ax2.get_xlim()
    ylim_ax2 = ax2.get_ylim()

    # Regions of the rest-frame lines
    xmins = np.array([3700, 4800, 6510])*(1+vals.z)
    xmaxs = np.array([3760, 5050, 6770])*(1+vals.z)
    for i in np.arange(len(xmins)):
        ax.axvspan(xmins[i], xmaxs[i], color='bisque', lw=0, zorder=0)

    # Connect the corners of the axvspan to the corners of the panel
    line1 = patches.ConnectionPatch((xmins[0], ylim_ax[1]),
            (xlim_ax2[0], ylim_ax2[0]), coordsA='data', coordsB='data', 
            axesA=ax, axesB=ax2, arrowstyle='-', color='bisque', lw=0.5)
    ax.add_artist(line1)

    line1 = patches.ConnectionPatch((xmaxs[0], ylim_ax[1]),
            (xlim_ax2[1], ylim_ax2[0]), coordsA='data', coordsB='data', 
            axesA=ax, axesB=ax2, arrowstyle='-', color='bisque', lw=0.5)
    ax.add_artist(line1)

    # Zoom-in of the middle lines
    ax2 = fig.add_subplot(gs[0, 1])
    panels(ax2, [4930-multi_width, 4930+multi_width], wl, flam, 0)
    panels(ax2, [4930-multi_width, 4930+multi_width], wl2, flam2, 1)
    plot_lines(ax2, 'oiii', vals.gc, lw=1)
    ax2.text(0.45, 0.95, '[O III]', ha='left', va='top', 
            transform=ax2.transAxes, color=vals.gc)
    plot_lines(ax2, 'hb', vals.rc, lw=0.3)
    ax2.text(0.15, 0.98, r'H$\beta$', ha='left', va='top', 
            transform=ax2.transAxes, color=vals.rc)
    ax2.set_xlabel("$\lambda_\mathrm{rest}$ $(\AA)$")

    # Get the corners of the panel
    xlim_ax2 = ax2.get_xlim()
    ylim_ax2 = ax2.get_ylim()

    # Connect the corners of the axvspan to the corners of the panel
    line1 = patches.ConnectionPatch((xmins[1], ylim_ax[1]),
            (xlim_ax2[0], ylim_ax2[0]), coordsA='data', coordsB='data', 
            axesA=ax, axesB=ax2, arrowstyle='-', color='bisque', lw=0.5)
    ax.add_artist(line1)

    line1 = patches.ConnectionPatch((xmaxs[1], ylim_ax[1]),
            (xlim_ax2[1], ylim_ax2[0]), coordsA='data', coordsB='data', 
            axesA=ax, axesB=ax2, arrowstyle='-', color='bisque', lw=0.5)
    ax.add_artist(line1)

    # Zoom-in of the right lines
    ax2 = fig.add_subplot(gs[0, 2])
    panels(ax2, [6640-multi_width, 6640+multi_width], wl, flam, 0)
    panels(ax2, [6640-multi_width, 6640+multi_width], wl2, flam2, 1)
    plot_lines(ax2, 'ha', vals.rc, lw=0.4)
    ax2.text(0.21, 0.98, r'H$\alpha$', ha='left', va='top', 
            transform=ax2.transAxes, color=vals.rc)
    plot_lines(ax2, 'nii', vals.gc, lw=1)
    ax2.text(0.3, 0.95, r'[N II]', ha='left', va='top', 
            transform=ax2.transAxes, color=vals.gc)
    plot_lines(ax2, 'sii', 'grey', lw=1)
    ax2.text(6725, 0.95, r'[S II]', ha='center', va='top', color='grey')
    ax2.set_xlabel("$\lambda_\mathrm{rest}$ $(\AA)$")

    # Get the corners of the panel
    xlim_ax2 = ax2.get_xlim()
    ylim_ax2 = ax2.get_ylim()

    # Connect the corners of the axvspan to the corners of the panel
    line1 = patches.ConnectionPatch((xmins[2], ylim_ax[1]),
            (xlim_ax2[0], ylim_ax2[0]), coordsA='data', coordsB='data', 
            axesA=ax, axesB=ax2, arrowstyle='-', color='bisque', lw=0.5)
    ax.add_artist(line1)

    line1 = patches.ConnectionPatch((xmaxs[2], ylim_ax[1]),
            (xlim_ax2[1], ylim_ax2[0]), coordsA='data', coordsB='data', 
            axesA=ax, axesB=ax2, arrowstyle='-', color='bisque', lw=0.5)
    ax.add_artist(line1)

    # Get data assuming z=0
    wl, flam, eflam = load_spec_1()
    wl2, flam2, eflam2 = load_spec_2()
    wl_lines = get_rest_wl()

    # Zoom-in of H-alpha at z=0 (6563 AA)
    ax2 = fig.add_subplot(gs[2, 0])
    panels(ax2, [6560-single_width, 6560+single_width], wl, flam, 0)
    panels(ax2, [6560-single_width, 6560+single_width], wl2, flam2, 1)
    plot_lines(ax2, 'ha', vals.rc)
    ax2.text(0.47, 0.95, r'H$\alpha$', ha='left', va='top', 
            transform=ax2.transAxes, color=vals.rc)
    #ax.set_xticks([6520, 6550, 6580])
    #ax.set_xticklabels([6520, 6550, 6580])
    ax2.set_xlabel("$\lambda_\mathrm{obs}$ ($\AA$)")

    # Get the corners of the panel
    xlim_ax2 = ax2.get_xlim()
    ylim_ax2 = ax2.get_ylim()

    # Connect the corners of the axvspan to the corners of the panel
    line1 = patches.ConnectionPatch((ls[0]-single_width, ylim_ax[0]),
            (xlim_ax2[0], ylim_ax2[1]), coordsA='data', coordsB='data', 
            axesA=ax, axesB=ax2, arrowstyle='-', color='lavender', lw=0.5)
    ax.add_artist(line1)

    line1 = patches.ConnectionPatch((ls[0]+single_width, ylim_ax[0]),
            (xlim_ax2[1], ylim_ax2[1]), coordsA='data', coordsB='data', 
            axesA=ax, axesB=ax2, arrowstyle='-', color='lavender', lw=0.5)
    ax.add_artist(line1)

    # Zoom-in of He I 5876
    ax2 = fig.add_subplot(gs[2, 1])
    panels(ax2, [5876-single_width, 5876+single_width], wl, flam, 0)
    panels(ax2, [5876-single_width, 5876+single_width], wl2, flam2, 1)
    plot_lines(ax2, 'hei', vals.rc)
    ax2.text(0.38, 0.98, r'He I 5876', ha='center', va='top', 
            transform=ax2.transAxes, color=vals.rc)
    ax2.set_xlabel("$\lambda_\mathrm{obs}$ ($\AA$)")
    #ax.set_xticks([5400, 5700, 6200])
    #ax.set_xticklabels([5400, 5700, 6200])

    # Get the corners of the panel
    xlim_ax2 = ax2.get_xlim()
    ylim_ax2 = ax2.get_ylim()

    # Connect the corners of the axvspan to the corners of the panel
    line1 = patches.ConnectionPatch((ls[1]-single_width, ylim_ax[0]),
            (xlim_ax2[0], ylim_ax2[1]), coordsA='data', coordsB='data', 
            axesA=ax, axesB=ax2, arrowstyle='-', color='lavender', lw=0.5)
    ax.add_artist(line1)

    line1 = patches.ConnectionPatch((ls[1]+single_width, ylim_ax[0]),
            (xlim_ax2[1], ylim_ax2[1]), coordsA='data', coordsB='data', 
            axesA=ax, axesB=ax2, arrowstyle='-', color='lavender', lw=0.5)
    ax.add_artist(line1)

    # Zoom-in of He II 4686
    ax2 = fig.add_subplot(gs[2, 2])
    panels(ax2, [4686-single_width, 4686+single_width], wl, flam, 0)
    panels(ax2, [4686-single_width, 4686+single_width], wl2, flam2, 1)
    plot_lines(ax2, 'heii', vals.rc)

    # Plot O II in the observer frame
    for l in wl_lines['oii']:
        ax2.axvline(l*(1+vals.z), lw=1, ymin=0, ymax=0.1, color=vals.gc)
    ax2.text(0.25, 0.01, r'[O II]', ha='left', va='bottom', 
            transform=ax2.transAxes, color=vals.gc)
    
    ax2.text(0.36, 0.98, r'He II 4686', ha='center', va='top', 
            transform=ax2.transAxes, color=vals.rc)
    ax2.set_xlabel("$\lambda_\mathrm{obs}$ ($\AA$)")
    #ax.set_xticks([4660, 4680, 4700, 4720])
    #ax.set_xticklabels([4660, 4680, 4700, 4720])

    # Get the corners of the panel
    xlim_ax2 = ax2.get_xlim()
    ylim_ax2 = ax2.get_ylim()

    # Connect the corners of the axvspan to the corners of the panel
    line1 = patches.ConnectionPatch((ls[2]-single_width, ylim_ax[0]),
            (xlim_ax2[0], ylim_ax2[1]), coordsA='data', coordsB='data', 
            axesA=ax, axesB=ax2, arrowstyle='-', color='lavender', lw=0.5)
    ax.add_artist(line1)

    line1 = patches.ConnectionPatch((ls[2]+single_width, ylim_ax[0]),
            (xlim_ax2[1], ylim_ax2[1]), coordsA='data', coordsB='data', 
            axesA=ax, axesB=ax2, arrowstyle='-', color='lavender', lw=0.5)
    ax.add_artist(line1)

    # Zoom-in of He II 4686
    plt.tight_layout()
    #plt.show()
    plt.savefig("spec.eps", dpi=300, bbox_inches='tight', pad_inches=0)
    plt.close()


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
    


