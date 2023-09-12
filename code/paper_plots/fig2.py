""" Combined Figure 2 

Required per Nature's editorial process """

import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "sans-serif"
from fig2_opt_lc import *
from fig2_mm_xray_comparison import *

# Initialize
figwidth_mm = 183 # Nature standard
figwidth_in = (figwidth_mm/10)/2.54 # in inches
fig = plt.figure(figsize=(figwidth_in, figwidth_in/1.2))

##  Plot full optical LC 

# Get data
dat = get_full_opt()

# Get the approximate time of explosion
t0_str = Time(vals.t0, format='jd').isot.replace("T", " ").split('.')[0]

# Get the time where you split the axes
tsplit = 20

# [left, bottom, width, height]
ax1 = fig.add_axes([0.1, 0.5, 0.3, 0.4])
ax2 = fig.add_axes([0.4, 0.5, 0.6, 0.4], sharey=ax1)

for ax in [ax1, ax2]:
    plot_22tsd(ax, show='apparent')
    plot_nonflare_epochs(ax, dat)

# Left panel: comparisons
ax = ax1
shift = np.abs(Planck18.distmod(z=0.0141).value-vals.dm)
plot_18cow(ax, show='apparent', offset=shift, single_col=True)
shift = np.abs(Planck18.distmod(z=0.2442).value-vals.dm)
plot_20xnd(ax, show='apparent', offset=shift)
shift = np.abs(Planck18.distmod(z=0.0085).value-vals.dm)
plot_98bw(ax, show='apparent', offset=shift)
shift = np.abs(Planck18.distmod(z=0.677).value-vals.dm)
#plot_sn2011kl(ax, show='apparent', offset=shift)
shift = np.abs(Planck18.distmod(z=0.1353).value-vals.dm)
plot_at2020mrf(ax, show='apparent', offset=shift-1.5)

# Left panel: spec epochs
plot_spec_epochs(ax)

# Left panel: formatting
ax.set_yticks([19,20,21,22,23])
ax.set_yticklabels([19,20,21,22,23])
ax.set_ylabel(r"$m_\mathrm{opt}$ (AB)", fontsize=10,
        fontname='sans-serif')
ax.set_xlim(-6,tsplit)

# Right panel: comparisons and epochs
ax = ax2
plot_spec_epochs(ax)
ax.text(16, 23.45, 'Opt. spec.', fontsize=8, ha='right', c='k')
plot_xray_epochs(ax)
ax.text(24.3, 23.25, 'X-ray', fontsize=8, ha='right', c='grey')
plot_radio_epochs(ax)
ax.text(26.7, 23, 'Radio', fontsize=8, ha='right', c='k')
plot_flare_epochs(ax, dat)
ax.text(26.6, 19.1, '$i$-band flare', fontsize=8, ha='right', c=vals.ic)
ax.text(25.8, 20.1, '$r$-band flare', fontsize=8, ha='right', c=vals.rc)
ax.text(68, 21.1, '$w$-band/clear flare', fontsize=8, ha='right', c=vals.wc)
ax.text(97, 20, '$g$-band flare', fontsize=8, ha='right', c=vals.gc)
shift = np.abs(Planck18.distmod(z=0.0085).value-vals.dm)
plot_98bw(ax, show='apparent', offset=shift)
shift = np.abs(Planck18.distmod(z=0.677).value-vals.dm)
shift = np.abs(Planck18.distmod(z=0.1353).value-vals.dm)
plot_at2020mrf(ax, show='apparent', offset=shift-1.5)

# Formatting of the right axis
ax.set_xlim(tsplit, 210)
ax.set_xscale('log')
ax.set_ylim(23.5, 18.8)
ax.set_xticks([30,40,50,70,100,200])
ax.set_xticklabels([30,40,50,70,100,200])

# Make a second y-axis
axtwin = ax.twinx()
axtwin.set_ylabel(
        "$M_\mathrm{opt}$ (AB)", fontsize=10, rotation=270, labelpad=15.0)
y_f = lambda y_i: y_i-vals.dm+2.5*np.log10(1+vals.z)
ymin, ymax = ax.get_ylim()
axtwin.set_ylim((y_f(ymin), y_f(ymax)))
axtwin.tick_params(axis='y', labelsize=10)
axtwin.set_yticks([-21,-20,-19,-18])
axtwin.set_yticklabels([-21,-20,-19,-18])
axtwin.plot([],[])

# For the legend
#ax.plot([100,100],[150,150],ls='-',c='grey',label='SN2011kl $g\'$')
ax.plot([100,100],[150,150],ls='--',c='k',label='AT2018cow $g$',lw=0.5)
#ax.plot([100,100],[150,150],ls='--',c=vals.rc,label='AT2018cow $r$',lw=0.5)
ax.plot([100,100],[150,150],ls='-',c='k',label='AT2020xnd $g$',lw=0.5)
ax.plot([100,100],[150,150],ls='-.',c='k',label='AT2020mrf $g$',lw=1)
#ax.plot([100,100],[150,150],ls='-',c=vals.rc,label='AT2020xnd $r$',lw=0.5)
ax.plot([100,100],[150,150],ls=':',c='k',label='SN1998bw $R_C$')
ax.scatter(0,0,marker='s',c=vals.gc,edgecolor='k',label='AT2022tsd $g$')
ax.scatter(0,0,marker='o',c=vals.rc,edgecolor='k',label='AT2022tsd $r$')
ax.scatter(0,0,marker='D',c=vals.ic,edgecolor='k',label='AT2022tsd $i$')
ax.scatter(0,0,marker='>',c='lightgrey',edgecolor='k',label='AT2022tsd $w$')
ax.scatter(0,0,marker='<',c=vals.oc,edgecolor='k',label='AT2022tsd $o$')
fig.legend(fontsize=8, bbox_to_anchor=(0.5,1.05), loc='upper center',
           ncol=5, handletextpad=0.1, columnspacing=0.8)

# Formatting of the whole fig
fig.text(0.5, 0, r'$\Delta t_\mathrm{obs}$ (days)', ha='center')
ax1.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
axtwin.spines['left'].set_visible(False)
ax1.tick_params(labelright='off')
ax2.yaxis.tick_right()
#axarr[1].set_yticks([])
axtwin.yaxis.tick_right()

fig.subplots_adjust(wspace=0)

# Millimeter panel
ax = fig.add_axes([0.1, 0.1, 0.4, 0.4])
run(ax)
ax.set_ylabel(
        r"$L_{\nu}$ ($\nu \gtrsim 100$ GHz; erg$\,$s$^{-1}\,$Hz$^{-1}$)",
        fontsize=10)
ax.tick_params(axis='both', labelsize=10)
ax.set_xlim(0.7, 300)
ax.set_ylim(1E25, 2E33)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r"$\Delta t_\mathrm{obs}$ (d)", fontsize=10)

# X-ray panel
ax = fig.add_axes([0.5, 0.1, 0.4, 0.4])
create_xray_panel(ax)

plt.show()
