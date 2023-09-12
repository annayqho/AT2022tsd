import matplotlib.pyplot as plt
import sys
sys.path.append("..")
import vals
import numpy as np
from matplotlib import rcParams
import matplotlib.patches as patches
rcParams['font.family'] = 'sans-serif'
from astropy.time import Time
from astropy.cosmology import Planck15
import pandas as pd
import cmasher as cmr
from ztfquery import marshal
from read_kann_lc import *

cowms=3

tde_col = vals.tde_col
lgrb_col = vals.lgrb_col
llgrb_col = vals.llgrb_col
iacol = 'grey'
slsncol = vals.sn_col
cccol = vals.sn_col
cowcol = vals.cow_col


def get_z(name):
    """ Get afterglow redshifts """
    dat = pd.read_csv("afterglows.txt")
    z = dat['z'][dat['Name']==name].values[0]
    return z



def calc_Lpeak(Mni, tpeak):
    """
    Peak luminosity as a function of nickel mass and rise time

    Mni: nickel mass in units of Msun
    """
    tauni = 8.8
    tauco = 113.6
    Lpeak = 2E43 * (Mni) * (3.9*np.exp(-tpeak/tauni) + \
          0.678*(np.exp(-tpeak/tauco)-np.exp(-tpeak/tauni)))
    return Lpeak


def calc_Mej(tpeak):
    """
    The ejecta mass corresponding to rise time (Eq 2 of Kasen et al. 2017)
    We're assuming kappa = 0.1, V9=1

    Mej: in units of Msun
    v9: in units of 1E9 cm/s
    tpeak: in units of days
    """
    v9 = 1.2 # 12,000 km/s
    Mej = ((tpeak*v9**(1/2))/14.5)**2
    return Mej


def plot_ztf(ax, background=False, shrink=1, text=True):
    """ Plot the ZTF sample """

    # Read table
    a = pd.read_csv("../../data/basic_info.csv")
    b = pd.read_csv("../../data/timescales.txt")

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

    i = 0
    for clname in np.unique(cl):
        if clname!='Unknown':
            choose = cl == clname

            # Scale by redshift
            x = dur[choose]/(1+z[choose])
            ex = edur[choose]/(1+z[choose])

            y = M[choose]
            ey = eM[choose]

            zorder = 10
            if clname=='?':
                label = 'unclassified'
                zorder=0
            if np.logical_or.reduce((clname=='II', clname=='IIb', clname=='Ic', 
                clname=='Ib', clname=='Ibn', clname=='Ic-BL', clname=='IIn/Ibn')):
                clname = 'II'
                label = 'Core-collapse SN'
            else:
                label=clname

            # Go onto plotting
            if background is False:
                if clname=='AT2018cow-like':
                    col = cowcol
                    m = 'D'
                    s = 10
                    label = 'LFBOT'
                else:
                    col = cccol
                    m = '>'
                    s = 10
                for j,name in enumerate(names[choose]):
                    ax.scatter(
                            x[j], y[j], label=None, c=col, 
                            marker=m, zorder=zorder, s=s)
                    if text:
                        if name=='ZTF18abcfcoo':
                            ax.text(x[j]*1.07, y[j], 'AT2018cow', fontsize=8,
                                    ha='left', va='bottom', c=cowcol)
                        if name=='ZTF18abvkwla':
                            ax.text(x[j]*1.07, y[j], 'AT2018lug', fontsize=8,
                                    ha='left', va='bottom', c=cowcol)
                        if name=='ZTF20acigmel':
                            ax.text(x[j]*1.07, y[j], 'AT2020xnd', fontsize=8,
                                    ha='left', va='bottom', c=cowcol)

    c = cowcol
    if background:
        c = 'lightgrey'
    # Plot label
    #if text:
    #    ax.text(5, -22.5, "LFBOT", va='center', ha='right', color=c)


def plot_comparison(ax):
    """ Plot the comparison sample from the literature """

    # Read table
    b = pd.read_table("../comparison_sample.txt", delimiter='&')
    names = b['Name'].values
    z = b['Redshift'].values
    cl = b['Class'].values
    
    peak = np.array([val.split('pm')[0] for val in b['Peak']]).astype(float)
    epeak = np.array([val.split('pm')[1] for val in b['Peak']]).astype(float)

    grise = b['Rise'].values.astype(str)
    gfade = b['Fade'].values.astype(str)

    egrise = np.copy(grise)
    egfade = np.copy(gfade)

    # For events with pm, can get error easily
    todo = np.array(['pm' in val for val in grise])
    grise[todo] = np.array([val.split('pm')[0] for val in grise[todo]])
    todo = np.array(['pm' in val for val in egrise])
    egrise[todo] = np.array([val.split('pm')[1] for val in egrise[todo]])
    todo = np.array(['pm' in val for val in gfade])
    gfade[todo] = np.array([val.split('pm')[0] for val in gfade[todo]])
    todo = np.array(['pm' in val for val in egfade])
    egfade[todo] = np.array([val.split('pm')[1] for val in egfade[todo]])

    # Munge the rises of the range events to have error bars
    tofix = np.array(['--' in val for val in grise])
    min_val = np.array([val.split('--')[0] for val in grise[tofix]]).astype(float)
    max_val = np.array([val.split('--')[1] for val in grise[tofix]]).astype(float)
    grise[tofix] = (min_val+max_val)/2
    egrise[tofix] = max_val-(min_val+max_val)/2

    tofix = np.array(['--' in val for val in gfade])
    min_val = np.array([val.split('--')[0] for val in gfade[tofix]]).astype(float)
    max_val = np.array([val.split('--')[1] for val in gfade[tofix]]).astype(float)
    gfade[tofix] = (min_val+max_val)/2
    egfade[tofix] = max_val-(min_val+max_val)/2

    # Now we're done with the fades
    gfade = gfade.astype(float)
    egfade = egfade.astype(float)

    # Munge the rises of the limit events
    tofix = np.array(['<' in val for val in grise])
    limval = np.array([val[1:] for val in grise[tofix]]).astype(float)
    grise[tofix] = limval
    egrise[tofix] = 0.0

    # Now we're done with the rises
    grise = grise.astype(float)
    egrise= egrise.astype(float)

    # Durations 
    dur = grise+gfade
    edur = np.sqrt(egrise**2+egfade**2)
   
    # Plot
    ax.errorbar(
            dur, peak, xerr=edur, yerr=epeak,
            label='Lit. Sample', c='grey', fmt='+', zorder=0)

    # Put arrows on events with limits
    islim = egrise==0
    for ii,xval in enumerate(grise[islim]):
        ax.arrow(xval, peak[islim][ii], -xval/10, 0, color='grey',
                head_length=xval/30,
                head_width=-peak[islim][ii]/100, length_includes_head=True, zorder=0)


def plot_afterglow(ax, name, x, y):
    """ Plot single event
    name: so that I can get the redshift
    x: duration in observer-frame
    y: luminosity in rest-frame g
    """
    z = float(get_z(name))
    xplt = x/(1+z)
    ax.scatter(x, y, c=lgrb_col, marker='P', s=20)
    if name!='ZTF19abvizsw':
        ax.arrow(x, y, -x/3, 0, length_includes_head=True,
                 head_width=np.abs(y/3), head_length=x/8, color=lgrb_col)
        ax.arrow(x, y, 0, np.abs(y)/1.3, length_includes_head=True,
                 head_width=np.abs(x/3.5), head_length=np.abs(y/2.5),
                 color=lgrb_col)


def plot_afterglows(ax):
    """ Add LGRB afterglows from the Kann sample """
    lc = load_kann_lc()
    for key,value in lc.items():
        z = float(value['z'])
        dm = Planck15.distmod(z=z).value
        if len(value['t'])>1:
            t = value['t'].astype(float) / (1+z)
            mag = value['mag'].astype(float)
            # Upon visual inspection, these have well-resolved peaks 
            # but observations only started after 100s
            keep_list = ['070419A', '070612A', '60707', '70802', '80710']
            if np.logical_or(sum(t<100/86400)>0, key in keep_list):
                #fig = plt.figure()
                # plt.plot(t, mag-dm)
                # plt.gca().invert_yaxis()
                # plt.xscale('log')
                # plt.axvline(x=100/86400, lw=0.5)
                # plt.title(key)
                # plt.savefig("afterglow_lc/%s.png" %key)
                # plt.close()
                peakm = np.min(mag)
                peakM = peakm-dm-2.5*np.log10(1+z)
                dur = min(t[mag>peakm+0.75])
                if dur>1E-1:
                    print(key)
                ax.scatter(dur, peakM, c=lgrb_col, marker='+', s=20)


def plot_AT2019pim(ax):
    """ Add AT2019pim to the plot """
    # we use the TESS LC fit from Dan's paper
    y = 1.8E45
    x = 0.36
    ax.scatter(x, y, marker='*', s=40, edgecolor=lgrb_col, facecolor='white')
    ax.text(x/1.3, y*2, 'Orphan', fontsize=7, c=lgrb_col, ha='right')
    ax.text(x/1.2, y, '(AT2019pim)', fontsize=7, c=lgrb_col, ha='right')


def plot_snls(ax):
    x = 3.81+8.60
    y = -20.26
    ey = 0.03
    ax.errorbar(x, y, yerr=ey, fmt='D', c=cowcol, ms=cowms)
    return x,y


def plot_11qr(ax):
    x = 2.9+8.7
    y = -19.56
    ey = 0.08
    ex = np.sqrt(0.1**2+0.4**2)
    ax.errorbar(x, y, yerr=ey, xerr=ex, fmt='o', c='k', zorder=50, lw=0.3)
    #ax.text(x/1.01, y, 'PS1-11qr', ha='right', va='bottom', zorder=50)


def plot_iptf15ul(ax):
    ec,fc,msize,shape = get_cc()
    x = 1.53+3.72
    ex = np.sqrt(0.05**2+0.08**2)
    y = -21.2
    ey = 0.3
    ax.errorbar(
            x, y, yerr=ey, xerr=ex, fmt=shape['Ibn'], 
            mfc=fc['Ibn'], c=fc['Ibn'], zorder=50, ms=msize['Ibn'], label='SN Ibn')
    #ax.text(x/1.01, y, 'iPTF15ul (Ibn?)', ha='right', va='bottom', zorder=50)


def plot_iptf16asu(ax):
    x = 1.14+10.62
    ex = np.sqrt(0.13**2+0.5**2)
    y = -20.3
    ey = 0.1
    ax.errorbar(
            x, y, yerr=ey, xerr=ex, ms=4,
            fmt='>', c='grey', zorder=50)


def plot_des14S2anq(ax):
    x = 10.44
    y = -16.19
    ey = 0.02
    ax.scatter(x, y, marker='o', c='k')
    ax.text(x*1.01, y*1.005, 'DES14S2anq', zorder=50)
    ax.arrow(x, y, -x/10, 0, color='k',
             head_length=x/30, head_width=-y/100, 
             length_includes_head=True)


def plot_des14S2plb(ax):
    x = 11.54
    y = -15.76
    ey = 0.13
    ax.scatter(x, y, marker='o', c='k')
    ax.text(x*1.01, y/1.005, 'DES14S2plb', zorder=50, ha='center', va='top')
    ax.arrow(x, y, -x/10, 0, color='k',
             head_length=x/30, head_width=-y/100, 
             length_includes_head=True)


def plot_des16S1dxu(ax):
    x = 11.16
    y = -16.04
    ey = 0.08
    ax.scatter(x, y, marker='o', c='k')
    ax.text(x*1.15, y*1.005, 'DES16S1dxu', zorder=50, ha='center')
    ax.arrow(x, y, -x/10, 0, color='k',
             head_length=x/30, head_width=-y/100, 
             length_includes_head=True)


def plot_sn1999cq(ax):
    """ Rise time is <3.9d, fade time is 10 days """
    x = 13.9
    y = -19.73
    ey = 0.1
    ec,fc,msize,shape = get_cc()
    ax.errorbar(x, y, yerr=ey, mfc=fc['Ibn'], mec=ec['Ibn'], marker=shape['Ibn'],
            ms=msize['Ibn'])
    ax.text(x, y/1.01, '99cq (Ibn)', zorder=50, ha='center', va='top')
    ax.arrow(x, y, -x/10, 0, color=fc['Ibn'],
             head_length=x/30, head_width=-y/100, 
             length_includes_head=True)


def plot_lsq12btw(ax):
    """ Rise time is <3.8d, fade time is 10 days """
    x = 13.9
    y = -19.3
    ey = 0.2
    ec,fc,msize,shape = get_cc()
    ax.errorbar(x, y, yerr=ey, mfc=fc['Ibn'], mec=ec['Ibn'], marker=shape['Ibn'],
            ms=msize['Ibn'])
    ax.text(x, y/1.01, '99cq (Ibn)', zorder=50, ha='center', va='top')
    ax.arrow(x, y, -x/10, 0, color=fc['Ibn'],
             head_length=x/30, head_width=-y/100, 
             length_includes_head=True)


def plot_sn2018kzr(ax):
    """ Rise time is <2d, fade time is 1.55pm0.23 days 

    So the minimum is 1.55-0.23=1.32 days, maximum is 1.55+0.23+2=3.78 days
    2.55 +/- 1.23
    """
    x = 2.55
    ex = 1.23
    y = -18.8
    ey = 0.08
    ec,fc,msize,shape = get_cc()
    ax.errorbar(x, y, xerr=ex, yerr=ey, mfc=fc['Ic'], mec=ec['Ic'], marker=shape['Ic'],
            ms=msize['Ic'])
    ax.text(x, y/1.01, '18kzr (Ic)', zorder=50, ha='center', va='top')


def plot_sn2019bkc(ax):
    """ Rise time is 5.28 pm 0.38, fade time is 2.22 pm 0.10 """
    x = 7.5
    ex = 0.39
    y = -17.16
    ey = 0.03
    ec,fc,msize,shape = get_cc()
    ax.errorbar(x, y, xerr=ex, yerr=ey, mfc=fc['Ic'], mec=ec['Ic'], marker=shape['Ic'], ms=msize['Ic'], c=fc['Ic'])
    ax.text(x, y, '19bkc (Ic)', zorder=50, ha='left', va='bottom')


def ibn_template(ax):
    """ Rise: -8 to -3... Fade: 4 to 10
    So, duration is between 7 and 18
    Peak is -18.9 to -19.8
    """
    ec,fc,msize,shape = get_cc()
    x = (7+18)/2
    ex = x-7
    y = -19.35
    ey = 0.55
    ax.errorbar(x, y, xerr=ex, yerr=ey, 
            mfc=fc['Ibn'], mec=ec['Ibn'], c=fc['Ibn'],
            marker=shape['Ibn'], ms=msize['Ibn'], zorder=50)
    ax.text(x, y/1.01, 'Ibn Template', zorder=50, ha='center', va='top')


def iib_first_peak(ax):
    """ 
    SN1993J and Fremling's object
    """
    ec,fc,msize,shape = get_cc()

    # ZTF18aalrxas
    x = 4.5
    y = -18.2
    ex = 0
    ey = 0
    ax.errorbar(x, y, xerr=ex, yerr=ey, 
            mfc=fc['IIb'], mec=ec['IIb'], c=fc['IIb'],
            marker=shape['IIb'], ms=msize['IIb'], zorder=50)
    ax.text(x, y/1.01, 'IIb, first peak', zorder=50, ha='center', va='top')

    # SN1993J 
    x = 3
    y = -17
    ex = 0
    ey = 0
    ax.errorbar(x, y, xerr=ex, yerr=ey, 
            mfc=fc['IIb'], mec=ec['IIb'], c=fc['IIb'],
            marker=shape['IIb'], ms=msize['IIb'], zorder=50)
    ax.text(x, y/1.01, 'IIb, first peak', zorder=50, ha='center', va='top')


def plot_bts(ax):
    """ Plot the BTS sample with consistent coloring """
    dat = pd.read_csv("../../data/bts.csv")
    dur = dat['duration'].values
    Mpk = dat['peakabs'].values
    cl = dat['type'].values
    names = dat['ZTFID'].values
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
    
    # Show them all in grey
    #ax.scatter(dur, Mpk, facecolor='grey', edgecolor='None', marker='s', alpha=0.3)

    choose = np.logical_or.reduce((cl=='SN II', cl=='SN IIb', cl=='SN Ic', cl=='SN Ib',
            cl=='SN Ibn', cl=='SN Ic-BL'))
    ax.scatter(
            dur[choose], Mpk[choose], 
            c=cccol, marker='>', zorder=2, s=10)
    #ax.text(17,-15.7,'Core-collapse',c=cccol, ha='right')#fontweight='bold',ha='right')
    #ax.text(13,-15.3,'SNe',c=cccol, ha='right')#, fontweight='bold'

    choose = cl=='SN Ia'
    ax.scatter(
            dur[choose], Mpk[choose], 
            c=iacol, marker='.', zorder=50, alpha=0.5, lw=0)

    choose = ['SLSN' in val for val in cl]
    ax.scatter(
            dur[choose], Mpk[choose],
            c=slsncol, marker='x', zorder=0, s=10)



def plot_panel(ax, zoom=False):
    """ Plot the luminosity-duration plot 
    If zoom=True, plot the zoom-in of the LFBOTs with the Mni line,
    labeling individual events.
    """

    # Plot BTS sources
    plot_bts(ax)

    # Plot ZTF sources
    if zoom:
        plot_ztf(ax, background=False, shrink=2, text=True)
    else:
        plot_ztf(ax, background=False, shrink=2, text=False)

    # Plot the first peaks of the LLGRB-SNe
    ax.scatter(1.0, -18.2, c=llgrb_col, marker='s', s=10)
    ax.scatter(1.4, -18.5, c=llgrb_col, marker='s', s=10)
    ax.scatter(1.4, -17.2, c=llgrb_col, marker='s', s=10)
    ax.scatter(1.2, -18.7, c=llgrb_col, marker='s', s=10)

    # Plot ZTF21abbxiyd
    x = 3
    y = -23.8
    ax.scatter(x, y, marker='D', s=10, c=cowcol)
    #ax.text(x/1.5,y/1.02,'AT2021ahuo', fontsize=labelsize,
    #        ha='center', va='top', color='red', fontweight='bold')

    # Plot CSS161010
    x = 5.5
    y = -20.7
    ax.errorbar(x, y, 
            label=None, mfc=cowcol, mec=cowcol,
            c=cowcol, fmt='D', ms=cowms, zorder=500)
    if zoom:
        ax.text(
                x/1.06, y, 'CSS161010', fontsize=7, 
                ha='right', va='center', c=cowcol)

    # Plot AT2020mrf
    ax.errorbar(7.1, -20, 
            label=None, mfc=cowcol, mec=cowcol,
            c=cowcol, fmt='D', ms=cowms, zorder=500)
    if zoom:
        ax.text(
                7.1/1.05, -20, 'AT2020mrf', fontsize=8, 
                ha='right', va='top', c=cowcol)

    # Plot GW170817
    if zoom==False:
        ax.scatter(0.6, -16, c='white', marker='*', 
                   edgecolor='grey', facecolor='white',s=40)
        ax.text(0.5, -16, 'AT2017gfo', ha='right', va='center', 
                c='grey', fontsize=8)

    # Plot DES16X1eho
    x = (1.28+2.53)/2 + 1.01
    y = -20.39
    ax.scatter(x, y, marker='D', c=cowcol, s=10)
    if zoom:
        ax.text(x*1.05, y, 'DES16X1eho', fontsize=8, 
                ha='left', va='bottom', c=cowcol)

    # Add individual objects from my Koala paper Table 1
    x,y = plot_snls(ax)
    if zoom:
        ax.text(x/1.05, y, 'SNLS04D4ec', fontsize=8, ha='right', va='top', c=cowcol)
    # Dougie: too luminous, don't show
    # SN 2011kl:
    x = 4.97+17.70
    ex = np.sqrt(1.20**2+5.82**2)
    y = -20.31
    ey = 0.13
    ax.errorbar(x, y, marker='*', xerr=ex, yerr=ey, ms=8,
                mec='grey', mfc='white', c='grey', lw=0.5)
    if zoom:
        ax.text(x/1.05, y, 'SN2011kl', fontsize=8, 
                ha='right', va='bottom', c='grey')

    # Plot AT2022tsd
    # 19.28 is the ext corr peak in g band
    # 19.74 is the value in r-band, which is g-band rest-frame
    Mpeak = 19.74-5*np.log10(vals.dL_mpc*1E6/10)+2.5*np.log10(1+vals.z)
    dur = 7.12
    edur = 2.57
    ax.errorbar(dur, Mpeak, xerr=edur, yerr=0.09, label=None, 
                c=cowcol, mec='k',
                fmt='D', ms=6, zorder=1000, lw=2)
    if zoom:
        ax.text(
                dur*1.03, Mpeak*1.001, 'AT2022tsd', va='bottom', ha='left', 
                color=cowcol, fontweight='bold')

    ax.set_ylim(-15,-37.5)
    if zoom:
        ax.set_ylim(-19.5, -21.7)
        ax.set_yticks([-20, -21])
        ax.set_yticklabels([-20, -21])
    
    # Luminosity axis
    ax2 = ax.twinx()
    # use the g-band effective wavelength: 4722.74 AA
    y_f = lambda y_i: 4E33*10**((y_i-4.77)/(-2.5))
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.plot([],[])
    ax2.set_yscale('log')
    ax2.tick_params(axis='both', labelsize=12)
    if zoom:
        ax2.set_ylabel(
                r"Peak $\nu L_\nu$ (erg s$^{-1}$)", fontsize=11, rotation=270, 
                labelpad=15.0)

    # Plot the AT2022tsd flares
    t = np.array([20,40,50,16,10,19,12,15,7,78,15,8])
    lum = np.array([1,2,2,6,4,5,0.8,2,2,3,0.8,0.6])
    ax2.scatter(t/(60*24*(1+vals.z)), lum*1E43, marker='*', 
                s=60, c=cowcol, edgecolor='k')
    if zoom==False:
        ax.text(1E-2, -21.6, 'AT2022tsd Flares', fontsize=10, c=cowcol,
                fontweight='bold', ha='center')

    # Plot GRB optical flash
    ax2.scatter(
            40/86400, 1E50, marker='*', edgecolor=lgrb_col, facecolor='white',
            zorder=500)
    if zoom==False:
        ax2.text(60/86400, 1E50, 'LGRB 080319B Prompt Flash', fontsize=8, c=lgrb_col)

    # Plot blazar flare
    ax2.scatter(30, 1E46, marker='*', edgecolor='grey', facecolor='white')

    # Plot AT2022cmc
    zcmc = 1.1930
    dm = Planck15.distmod(z=zcmc).value
    Mpeak = 19.25-dm+2.5*np.log10(1+zcmc)
    x = 0.31+0.89
    ax.scatter(x,Mpeak,marker='o',s=20, zorder=500,c=vals.tde_col)

    # Plot ASASSN-14ko 
    #ax2.scatter(10, 1E43, marker='*', edgecolor='grey', facecolor='white')
    
    if zoom==False:
        ax2.text(160, 3.0E46, 'S5 1803+784', fontsize=7, c='grey', ha='right')
        ax2.text(160, 1.5E46, 'Blazar Flare', fontsize=7, c='grey', ha='right')
        ax.text(x/1.1,Mpeak*1.01,'AT2022cmc',c=vals.tde_col,
                ha='left',va='bottom',fontsize=7)
        #ax2.text(160, 2.8E46, 'ASASSN-14ko', fontsize=8, c='grey', ha='right')
        #ax2.text(50, 1.5E43, 'Partial SMBH TDE?', fontsize=8, 
        #         c='grey', ha='right')

    if zoom==False:
        # Plot afterglows
        plot_afterglows(ax)
        plot_AT2019pim(ax2)
        ax.set_ylabel("$M_{g,\mathrm{peak}}$", fontsize=11)

    ax.set_xlim(7E-5,200)
    ax.set_xscale('log')
    if zoom:
        ax.set_xlim(2.5,40)
        ax.set_xticks([4,10,20,40])
        ax.set_xticklabels([4,10,20,40])
    ax.set_xlabel(r"$t_{1/2, \mathrm{optical}}$ (rest-frame days)", fontsize=11)
    ax.tick_params(axis='both', labelsize=12)

    if zoom:
        # Plot the Mej = Mni line
        # Plot the Mni=Mej limit
        tpeak = np.linspace(1, 300)
        Mej = calc_Mej(tpeak)
        Lpeak = calc_Lpeak(Mej, tpeak)
        ax2.plot(tpeak, Lpeak, c='k', ls='--')
        ax.text(11.4, -19.8, r'$M_{\mathrm{Ni}}>M_{\mathrm{ej}}$', rotation=60)
        ax.text(13.3, -19.7, r'$M_{\mathrm{Ni}}<M_{\mathrm{ej}}$', rotation=60)

        ax.scatter(0, 0, marker='>', c=cccol, label='CC SN')
        ax.scatter(0, 0, marker='P', c=lgrb_col, label='LGRB Afterglows')
        ax.scatter(0, 0, marker='x', c=slsncol, label='SLSN')
        ax.scatter(0, 0, marker='o', c=iacol, label='SN Ia')
        ax.scatter(0, 0, marker='D', c=cowcol, label='LFBOT')
        ax.scatter(0, 0, marker='s', c=llgrb_col, label='LLGRB-SN First Peak')
        ax.scatter(0, 0, marker='o', c=vals.tde_col, label='TDE')
