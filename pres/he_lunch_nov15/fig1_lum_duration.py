"""
Figure 1 in the paper.
Plot the duration vs. luminosity of various optical transients,
as measured by ZTF.
"""

# Import packages
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from astropy.cosmology import Planck15
import pandas as pd
from ztfquery import marshal
import cmasher as cm

# Import my own functions
import sys
from helpers import get_z

# Global variables
labelsize = 8
cols = cm.take_cmap_colors(
        'cmr.rainforest', 5, cmap_range=(0.0, 0.85), return_fmt='hex')
iacol = cols[4]
cccol = cols[1]
aftcol = cols[2]
cowcol = cols[3]
slsncol = cols[0]

# Functions

def plot_text(ax):
    """ Plot the text labels """
    ax.text(
            30,-23,'SLSNe',c=slsncol, fontweight='bold', fontsize=labelsize)
    ax.text(
            10,-20.3,'SN Ia',
            c=iacol, rotation=30, fontweight='bold', fontsize=labelsize)
    ax.text(
            3,-16,'Core-collapse SNe',
            c=cccol, ha='right',fontsize=labelsize,fontweight='bold')
    ax.text(
            2,-21,'LFBOT',
            c=cowcol, ha='right',fontsize=labelsize,fontweight='bold')
    ax.text(0.4, -19, 'Rel. SBO', c='grey', fontweight='bold',
            fontsize=labelsize)
    ax.text(
            0.1, -27, 'Orphan Afterglows', 
            color=aftcol, fontsize=labelsize, fontweight='bold')


def plot_cc(ax,x,y):
    """ Plot the CC SNe """
    ax.scatter(x, y, c=cccol, marker='o', zorder=2, s=5)


def plot_cow(ax,x,y,ztf=True):
    """ Plot the Cows """
    if ztf:
        ax.scatter(x, y, c=cowcol, marker='D', zorder=2, s=15,
                edgecolor='k', facecolor=cowcol)
    else:
        ax.scatter(x, y, c=cowcol, marker='D', zorder=2, s=15,
                edgecolor=cowcol, facecolor=cowcol)


def plot_afterglow(ax, name, x, y):
    """ Plot single event 
    name: so that I can get the redshift
    x: duration in observer-frame
    y: luminosity in rest-frame g
    """
    msize = 30
    z = float(get_z(name))
    xplt = x/(1+z)
    if name=='iPTF14yb':
        ax.scatter(x, y, edgecolor=aftcol, facecolor='white', marker='o', 
                s=msize, zorder=50)
        ax.text(
                x/1.2, y/1.2, 'iPTF14yb', c=aftcol, ha='right', 
                va='top', fontsize=8) 
    else:
        ax.scatter(x, y, c=aftcol, marker='o', s=msize, zorder=2,
                edgecolor='k', facecolor=aftcol)
    if name!='ZTF19abvizsw':
        ax.arrow(x, y, -x/2, 0, length_includes_head=True,
                 head_width=np.abs(y/3), head_length=x/6, color=aftcol) 
        ax.arrow(x, y, 0, np.abs(y)/1.3, length_includes_head=True,
                 head_width=np.abs(x/3.5), head_length=np.abs(y/2.5), 
                 color=aftcol) 


def plot_ztf(ax, background=False, shrink=1, text=True):
    """ Plot Rapidly Evolving Transients from ZTF (Ho+2021) """

    # Read table
    a = pd.read_csv("basic_info.csv")
    b = pd.read_csv("timescales.txt")

    # Get basic info
    names = a['Name'].values
    z = a['Redshift'].values
    cl = a['Class'].values
    m = a['mgpeak'].values
    em = a['emgpeak'].values

    # Get timescales
    name_key = b['name'].values
    grise = np.array([b['grise'].values[name_key==val][0] for val in names])
    gfade = np.array([b['gfade'].values[name_key==val][0] for val in names])

    # Munge the rises of the range events to have error bars
    tofix = np.array(['-' in val for val in grise])
    min_val = np.array([val.split('-')[0] for val in grise[tofix]]).astype(float)
    max_val = np.array([val.split('-')[1] for val in grise[tofix]]).astype(float)
    avgval = np.round((min_val+max_val)/2,2)
    errval = avgval-min_val
    repval = np.array(['%spm%s' %(i,j) for i,j in zip(avgval,errval)])
    grise[tofix] = repval

    # Munge the rises of the limit events
    tofix = np.array(['<' in val for val in grise])
    limval = np.array([val[1:] for val in grise[tofix]]).astype(float)
    grise[tofix] = np.array(['%spm%s' %(val,0) for val in limval])

    # Munge the fades of the limit events
    tofix = np.array(['>' in val for val in gfade])
    limval = np.array([val[1:] for val in gfade[tofix]]).astype(float)
    gfade[tofix] = np.array(['%spm%s' %(val,0) for val in limval])

    # Only keep events with known redshift
    keep = ~np.isnan(z)
    names = names[keep]
    z = z[keep].astype(float)
    cl = cl[keep]
    trise = np.array([float(val.split('pm')[0]) for val in grise[keep]])
    etrise = np.array([float(val.split('pm')[1]) for val in grise[keep]])
    tfade = np.array([float(val.split('pm')[0]) for val in gfade[keep]])
    etfade = np.array([float(val.split('pm')[1]) for val in gfade[keep]])

    M = m[keep].astype(float)-Planck15.distmod(z=z).value
    eM = em[keep].astype(float)

    # Durations
    dur = trise+tfade
    edur = np.sqrt(etrise**2+etfade**2)

    # Scale by redshift
    dur = dur/(1+z)
    edur = edur/(1+z)

    # Plot the CC SNe
    choose = np.logical_or.reduce((cl=='II', cl=='IIb', cl=='Ic', cl=='Ib',
            cl=='Ibn', cl=='Ic-BL', cl=='IIn/Ibn'))
    y = M[choose]
    plot_cc(ax, dur[choose], M[choose])

    # Plot the Cows
    choose = names=='ZTF18abcfcoo'
    plot_cow(ax, dur[choose], M[choose], ztf=False)
    choose = names=='ZTF20acigmel'
    plot_cow(ax, dur[choose], M[choose], ztf=True)
    choose = names=='ZTF18abvkwla'
    plot_cow(ax, dur[choose], M[choose], ztf=True)


def plot_bts(ax):
    """ Plot the BTS sample with consistent coloring """

    # Load the BTS measurements
    dat = pd.read_csv("bts.csv")
    dur = dat['duration'].values
    Mpk = dat['peakabs'].values
    cl = dat['type'].values
    names = dat['ZTFID'].values

    # Get rid of null values, limits
    keep = np.array([val!='-' for val in Mpk])
    dur = dur[keep]
    cl = cl[keep]
    names = names[keep]
    Mpk = Mpk[keep].astype(float)
    keep = np.array(['>' not in val for val in dur])
    dur = dur[keep].astype(float)
    Mpk = Mpk[keep]
    cl = cl[keep]
    names = names[keep]
    
    # Plot the core-collapse supernovae
    choose = np.logical_or.reduce(
            (cl=='SN II', cl=='SN IIb', cl=='SN Ic', 
             cl=='SN Ib', cl=='SN Ibn', cl=='SN Ic-BL'))
    ax.scatter(
            dur[choose], Mpk[choose], 
            c=cccol, marker='o', zorder=2, s=5)

    # Plot the Typa Ia supernovae
    choose = cl=='SN Ia'
    ax.scatter(
            dur[choose], Mpk[choose], 
            c=iacol, marker='.', zorder=0)

    # Plot the SLSNe
    choose = ['SLSN' in val for val in cl]
    ax.scatter(
            dur[choose], Mpk[choose], 
            c=slsncol, marker='x', zorder=0, s=10)


def plot_afterglows(ax):
    """ Add afterglows to the plot """

    # ZTF19abvizsw: we use the TESS LC fit from Dan's paper
    y = 1.8E45
    x = 0.36
    plot_afterglow(ax, 'ZTF19abvizsw', x, y)

    # ZTF20aajnksq: we use the first r-band detection
    # and for the duration, the time from estimated t0 to first det (2.4 hr)
    # + the approximate time to half-max assuming the power law (8.6hr)
    y = 9.4E45
    x = (11/24)
    plot_afterglow(ax, 'ZTF20aajnksq', x, y)

    # ZTF21aaeyldq: we use the first r-band detection
    # and for the duration, the time from estimated t0 to first det (14 min)
    # plus the time from peak to half-max using the power law (36 min)
    # so total time is 50 min = 50/60/24
    y = 3.8E46
    x = (50/60/24)
    plot_afterglow(ax, 'ZTF21aaeyldq', x, y)

    # ZTF21aayokph:
    # time from t0 to first det: 0.98d
    # estimated time of fade: 1.3d...2.25d
    y = 2.8E45
    x = 2.25
    plot_afterglow(ax, 'ZTF21aayokph', x, y)

    # iPTF14yb
    # time from GRB to first det: 14.7 minutes
    # time to from first det to half: .04d
    #y = 5.1E45
    #x = ((14.7/60/24) + 0.04)
    #plot_afterglow(ax, 'iPTF14yb', x, y)

    # AT2020kym
    # time from GRB to first det: 1.8hr
    # time from first det to half-max: 0.034d
    # y = 1.3E46
    # x = ((1.8/24) + 0.034)
    # plot_afterglow(ax, 'ZTF20abbiixp', x, y)

    # ZTF20acozryr
    # time from GRB to first det: 0.6d
    # time to half-max: 1 day
    # y = 1.6E45
    # x = 1.6
    # plot_afterglow(ax, 'ZTF20acozryr', x, y)

    # ZTF21aagwbjr
    # time from GRB to first det: 43 minutes
    # time to half-max: 137 minutes (from t0)
    # y = 6.4E45
    # x = (137-43)/60/24
    # plot_afterglow(ax, 'ZTF21aagwbjr', x, y)

    # ZTF21abfmpwn
    # time from GRB to first det: 9.7 hr
    # time to half-max: 15 hours (from t0)
    # y = 4.7E45
    # x = (15-9.7)/24
    # plot_afterglow(ax, 'ZTF21abfmpwn', x, y)



if __name__=="__main__":
    # Initialize the figure
    fig,ax = plt.subplots(1,1,figsize=(4,4))

    # Plot BTS sources
    plot_bts(ax)

    # Plot sources from my RET paper
    plot_ztf(ax, background=False, shrink=2, text=True)

    # Plot AT2020mrf
    plot_cow(ax, 7.1, -20)

    # Plot CSS161010
    plot_cow(ax, 5.5, -21.5, ztf=False)

    # Plot GW170817 
    #ax.scatter(0.6, -16, c='white', marker='s', edgecolor='k', facecolor='white')
    #ax.text(0.5, -16, 'AT2017gfo', ha='right', va='center', c='k', fontsize=9)

    # Add the Rel. SBO
    ax.scatter(1.0, -18.2, c='grey', marker='s', edgecolor='k', s=10)
    ax.scatter(1.4, -18.5, c='grey', marker='s', s=10)
    ax.scatter(1.4, -17.2, c='grey', marker='s', s=10)
    ax.scatter(1.2, -18.7, c='grey', marker='s', s=10)

    # Plot AT2022tsd g-band LC
    z = 0.256
    Mpeak = -20.2
    x = 9
    ax.scatter(x,Mpeak,facecolor=cowcol,marker='*',s=150, zorder=500,
            edgecolor='k')
    ax.text(x,Mpeak/1.03,'AT2022tsd', fontsize=labelsize, ha='center', 
            va='top', color=cowcol, fontweight='bold')

    ax.set_ylabel("$M_{g,\mathrm{peak}}$ (rest-frame)", fontsize=11)
    ax.set_ylim(-15,-28.5)

    # Twin axis for luminosity
    ax2 = ax.twinx()
    # use the g-band effective wavelength: 4722.74 AA
    y_f = lambda y_i: 4E33*10**((y_i-4.77)/(-2.5))
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.plot([],[])
    ax2.set_yscale('log')
    ax2.tick_params(axis='both', labelsize=10)
    ax2.set_ylabel(
            r"Peak $\nu L_\nu$ (erg s$^{-1}$)", fontsize=11, rotation=270, 
            labelpad=15.0)

    ax.set_zorder(ax2.get_zorder()+1)
    ax.patch.set_visible(False)

    # Plot the afterglows
    plot_afterglows(ax2)

    # Plot the text
    plot_text(ax)

    # Plot the lines with our search criteria
    #ax.axhline(y=-19.8,c='grey',ls='--',lw=0.5)
    #ax.axvline(x=8,c='grey',ls='--',lw=0.5)

    ax.set_xlim(0.01,200)
    ax.set_xscale('log')
    ax.set_xlabel("Rest-frame days above half-max", fontsize=11)
    ax.tick_params(axis='both', labelsize=10)

    #plt.tight_layout()
    plt.show()
    #plt.savefig("lum_time_optical.png", dpi=200, bbox_inches='tight', pad_inches=0.1)
    #plt.close()
