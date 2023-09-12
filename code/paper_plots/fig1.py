""" Combine the two panels of Fig 1 into a single figure

Final part of Nature editorial process...
"""

import matplotlib.gridspec as gridspec
from fig1_host_galaxy import *
from fig1_lum_duration import *

figwidth_mm = 183 # Nature standard
figwidth_in = (figwidth_mm/10)/2.54 # in inches
fig,axs = plt.subplots(2,2,figsize=(figwidth_in,figwidth_in))

# Luminosity-duration panels
axarr = axs[0,:]
ax = axarr[0]
plot_panel(ax)
ax = axarr[1]
plot_panel(ax, zoom=True)

# Get the corners of the second panel
xlim_ax2 = axarr[1].get_xlim()
ylim_ax2 = axarr[1].get_ylim()
width_ax2 = abs(xlim_ax2[1]-xlim_ax2[0])
height_ax2 = ylim_ax2[1]-ylim_ax2[0]

# Draw a rectangle in the first panel
rect = patches.Rectangle(
        (xlim_ax2[0], ylim_ax2[0]), width_ax2, height_ax2,
        linewidth=0.5, edgecolor='k', facecolor='none', zorder=100)
axarr[0].add_patch(rect)

# Connect the corners of the rectangle to the corners of the second panel
# Upper right to top left
line1 = patches.ConnectionPatch(
        (xlim_ax2[0]+width_ax2, ylim_ax2[0]+height_ax2),
        (xlim_ax2[0], ylim_ax2[0]+height_ax2),
        coordsA='data', coordsB='data', axesA=axarr[0], axesB=axarr[1],
        arrowstyle='-', color='grey', lw=0.5)
axarr[0].add_artist(line1)

# Lower right to bottom left
line1 = patches.ConnectionPatch(
        (xlim_ax2[0]+width_ax2, ylim_ax2[0]),
        (xlim_ax2[0], ylim_ax2[0]),
        coordsA='data', coordsB='data', axesA=axarr[0], axesB=axarr[1],
        arrowstyle='-', color='grey', lw=0.5)
axarr[0].add_artist(line1)

fig.legend(loc='upper center', fontsize=10, handletextpad=0.1,
      bbox_to_anchor=(0.52, 1.00),
      ncol=7, fancybox=True, columnspacing=0.4)
