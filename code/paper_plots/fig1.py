""" Combine the two panels of Fig 1 into a single figure

Final part of Nature editorial process...
"""

from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
import matplotlib.gridspec as gridspec
from fig1_host_galaxy import *
from fig1_lum_duration import *

figwidth_mm = 183 # Nature standard
figwidth_in = (figwidth_mm/10)/2.54 # in inches
fig,axs = plt.subplots(2,2,figsize=(figwidth_in,figwidth_in/1.2))

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

fig.legend(loc='upper center', fontsize=10, handletextpad=0.05,
      bbox_to_anchor=(0.51, 1.00), 
      ncol=7, fancybox=True, columnspacing=0.4)

# Plot the host galaxy
ax = axs[1,0]

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

ax.text(0.01, 0.99,"Keck/LRIS $ugI$",fontsize=12,transform=ax.transAxes,
    horizontalalignment='left', va='top', color='white')

# Mark image scale : 0.25 arcsec per pixel
x = 7
y = 10
x2 = x + 5/0.25
ax.plot((x,x2), (y,y), color=markcol, lw=2)
ax.text((x2+x)/2, y*1.1, "5''", color=markcol, fontsize=12,
        verticalalignment='bottom', horizontalalignment='center')
ax.text((x2+x)/2, y/1.1, "(21 kpc)", color=markcol, fontsize=12,
        verticalalignment='top', horizontalalignment='center')


# Turn off final panel
axs[1,1].set_visible(False)

axs[0,0].text(1.15, 1.15, 'a', transform=ax.transAxes,
      fontsize=11, fontweight='bold', va='top', ha='right')
axs[1,0].text(1.05, 1.0, 'b', transform=ax.transAxes,
      fontsize=11, fontweight='bold', va='top', ha='right')

# Save
#fig.subplots_adjust(wspace=0.5, hspace=0.5)
plt.tight_layout()
fig.subplots_adjust(top=0.93, wspace=0.4)
plt.savefig("fig1.eps", dpi=300, bbox_inches='tight', pad_inches=0)
plt.close()
#plt.show()
