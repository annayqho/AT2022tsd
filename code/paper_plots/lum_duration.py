import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
import vals
import numpy as np
from astropy.time import Time
from astropy.cosmology import Planck15
import pandas as pd
import cmasher as cmr
from ztfquery import marshal

iacol = '#003f5c'
#cccol = '#58508d'
aftcol = '#bc5090'
llgrbcol = '#ff6361'
cowcol = '#ffa600'

cols = cmr.take_cmap_colors(
    'cmr.rainforest', 5, cmap_range=(0.1, 0.9), return_fmt='hex')[::-1]

tde_col = cols[4]
lgrb_col = cols[0]
llgrb_col = cols[3]
sn_col = 'grey'
cow_col = cols[2]

cccol = sn_col
cowcol = cow_col


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
    v9 = 1.5 # 15,000 km/s
    Mej = ((tpeak*v9**(1/2))/14.5)**2
    return Mej


def plot_ztf(ax, background=False, shrink=1, text=True):
    """ Plot the ZTF sample """

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
                    s = 5
                else:
                    col = cccol
                    m = '>'
                    s = 3
                for j,name in enumerate(names[choose]):
                    print(name)
                    ax.errorbar(
                            x[j], y[j], xerr=ex[j], yerr=ey[j], 
                            label=None, mfc=col, mec=col,
                            c=col, fmt=m, zorder=zorder, ms=s)

    c = cowcol
    if background:
        c = 'lightgrey'
    # Plot label
    if text:
        ax.text(5, -22.5, "LFBOT", va='center', ha='right', color=c)


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


def plot_snls(ax):
    x = 3.81+8.60
    y = -20.26
    ey = 0.03
    ax.errorbar(x, y, yerr=ey, fmt='D', c=cowcol)#, label='SNLS unclassified')


def plot_10ah(ax):
    x = 1+6.3
    y = -17.59
    ey = 0.11
    ex = np.sqrt(0.1**2+0.6**2)
    ax.errorbar(x, y, yerr=ey, xerr=ex, fmt='o', c='k', zorder=50, lw=0.3)
    #ax.text(x*1.01, y/1.005, 'PS1-10ah', ha='left', va='top', zorder=50)


def plot_10bjp(ax):
    x = 1+7.7
    y = -18.34
    ey = 0.11
    ex = np.sqrt(0.1**2+0.6**2)
    ax.errorbar(x, y, yerr=ey, xerr=ex, fmt='o', c='k', zorder=50, lw=0.3)
    #ax.text(x/1.01, y, 'PS1-10bjp', ha='right', va='top', zorder=50)


def plot_11qr(ax):
    x = 2.9+8.7
    y = -19.56
    ey = 0.08
    ex = np.sqrt(0.1**2+0.4**2)
    ax.errorbar(x, y, yerr=ey, xerr=ex, fmt='o', c='k', zorder=50, lw=0.3)
    #ax.text(x/1.01, y, 'PS1-11qr', ha='right', va='bottom', zorder=50)


def plot_12bb(ax):
    x = (6.3+1.8+6.3)/2
    ex = x-6.3
    y = -16.97
    ey = 0.12
    ax.errorbar(x, y, yerr=ey, xerr=ex, fmt='o', c='k', zorder=50, label='PS1 unclassified', lw=0.3)
    #ax.text(x/1.01, y, 'PS1-12bb', ha='right', va='bottom', zorder=50)


def plot_12bv(ax):
    """ <2.2; 3-9; minimum: 3...maximum: 2.2+9=11.2 """
    x = (3+11.2)/2
    ex = x-3
    y = -19.49
    ey = 0.07
    ax.errorbar(x, y, yerr=ey, xerr=ex, fmt='o', c='k', zorder=50, lw=0.3)
    #ax.text(x/1.01, y, 'PS1-12bv', ha='right', va='bottom', zorder=50)


def plot_12brf(ax):
    """ <1; 8.8+/-0.6 ... so from 8.2 to 10.4 """
    x = (8.2+10.4)/2
    ex = x-8.2
    y = -18.43
    ey = 0.08
    ax.errorbar(x, y, yerr=ey, xerr=ex, fmt='o', c='k', zorder=50, lw=0.3)
    #ax.text(x*1.01, y, 'PS1-12brf', ha='left', va='bottom', zorder=50)


def plot_ptf12ldy(ax):
    ec,fc,msize,shape = get_cc()
    x = 3.34+7.57
    ex = np.sqrt(0.17**2+0.29**2)
    y = -19.20
    ey = 0.02
    ax.errorbar(
            x, y, yerr=ey, xerr=ex, fmt=shape['Ibn'], 
            mfc=fc['Ibn'], c=fc['Ibn'], zorder=50, ms=msize['Ibn'], lw=0.3)
    #ax.text(x/1.05, y/1.005, 'PTF12ldy (Ibn)', ha='left', va='top', zorder=50)


def plot_13dwm(ax):
    """ <3.0; 3-7. So range is 3-10 """
    x = (3+10)/2
    ex = x-3
    y = -17.63
    ey = 0.13
    ax.errorbar(x, y, yerr=ey, xerr=ex, fmt='o', c='k', zorder=50, lw=0.3)
    #ax.text(x/1.01, y, 'PS1-13dwm', ha='right', va='bottom', zorder=50)


def plot_lsq13ccw(ax):
    ec,fc,msize,shape = get_cc()
    x = 1.39+3.86
    ex = np.sqrt(0.10**2+0.31**2)
    y = -18.4
    ey = 0.2
    ax.errorbar(
            x, y, yerr=ey, xerr=ex, fmt=shape['IIn/Ibn?'], 
            mfc=fc['IIn/Ibn?'], c=fc['IIn/Ibn?'], 
            zorder=50, ms=msize['IIn/Ibn?'])
    ax.text(x/1.01, y, 'LSQ13ccw (IIn/Ibn?)', 
            ha='right', va='bottom', zorder=50)


def plot_iptf14aki(ax):
    ec,fc,msize,shape = get_cc()
    x = 3.34+7.58
    ex = np.sqrt(0.17**2+0.30**2)
    y = -19.30
    ey = 0.03
    ax.errorbar(
            x, y, yerr=ey, xerr=ex, fmt=shape['Ibn'], 
            mfc=fc['Ibn'], c=fc['Ibn'], zorder=50, ms=msize['Ibn'])
    #ax.text(x/1.01, y, 'iPTF14aki (Ibn)', ha='right', va='bottom', zorder=50)


def plot_iptf15akq(ax):
    ec,fc,msize,shape = get_cc()
    x = 3.13+8.86
    ex = np.sqrt(0.61**2+0.80**2)
    y = -18.62
    ey = 0.31
    ax.errorbar(
            x, y, yerr=ey, xerr=ex, fmt=shape['Ibn'], 
            mfc=fc['Ibn'], c=fc['Ibn'], zorder=50, ms=msize['Ibn'])
    #ax.text(x/1.01, y, 'iPTF15akq (Ibn)', ha='right', va='bottom', zorder=50)


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


def plot_ksn(ax):
    x = 1.3+5.9
    y = -18.8
    # need to set
    ex = 0
    ey = 0
    ax.errorbar(x, y, yerr=ey, xerr=ex, fmt='*', c='k', zorder=50, ms=10)
    ax.text(x/1.01, y, 'KSN2015K', ha='right', va='top', zorder=50)


def plot_des(ax):
    """ Plot all DES objects from file """
    dat = pd.read_csv("../des_timescales.txt")
    names = dat['Name'].values
    grise = dat['Rise'].values
    gfade = dat['Fade'].values
    Mraw = dat['Peak'].values

    keep = np.loadtxt("../keep_des.txt", dtype=str)
    tokeep = np.array([val in keep for val in names])
    names = names[tokeep]
    grise =grise[tokeep]
    gfade = gfade[tokeep]
    Mraw = Mraw[tokeep]

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

    # Munge the peak mags
    M = np.array([val.split('+/-')[0] for val in Mraw]).astype(float)
    eM = np.array([val.split('+/-')[1] for val in Mraw]).astype(float)

    # Durations
    trise = np.array([val.split('pm')[0] for val in grise]).astype(float)
    etrise = np.array([val.split('pm')[1] for val in grise]).astype(float)
    tfade = np.array([val.split('pm')[0] for val in gfade]).astype(float)
    etfade = np.array([val.split('pm')[1] for val in gfade]).astype(float)
    dur = trise+tfade
    edur = np.sqrt(etrise**2+etfade**2)

    ax.errorbar(dur, M, yerr=eM, xerr=edur, fmt='s', c='k', 
            zorder=10, label='DES unclassified', lw=0.3, ms=5)



def plot_iptf16asu(ax):
    x = 1.14+10.62
    ex = np.sqrt(0.13**2+0.5**2)
    y = -20.3
    ey = 0.1
    ax.errorbar(
            x, y, yerr=ey, xerr=ex, ms=msize['Ic-BL'],
            fmt=shape['Ic-BL'], c=fc['Ic-BL'], zorder=50)
    ax.text(x*1.01, y*1.005, 'iPTF16asu (Ic-BL)', 
            ha='right', va='bottom', zorder=50, c='cornflowerblue')


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
    dat = pd.read_csv("bts.csv")
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
    ax.text(17,-15.7,'Core-collapse',c=cccol, ha='right')#fontweight='bold',ha='right')
    ax.text(13,-15.3,'SNe',c=cccol, ha='right')#, fontweight='bold'

    choose = cl=='SN Ia'
    ax.scatter(
            dur[choose], Mpk[choose], 
            c=iacol, marker='.', zorder=0, alpha=0.5, lw=0)
    ax.text(14.8,-19.6,'SN Ia',c=iacol, rotation=35)#, fontweight='bold')

    choose = ['SLSN' in val for val in cl]
    ax.scatter(
            dur[choose], Mpk[choose], 
            c='grey', marker='x', zorder=0, s=10)
    ax.text(40,-22.5,'SLSN',c='grey')#, fontweight='bold')



def plot_afterglows(ax):
    """ Add afterglows to the plot """
    x = 0.3
    y = -24.4
    ax.scatter(x, y, c=aftcol, edgecolor='k') 
    #ax.text(0.3, -24.4, 'ZTF19abvizsw', c=aftcol) 

    x = 0.64
    y = -26
    ax.scatter(x, y, c=aftcol, zorder=10) 
    ax.arrow(x, y, -x/3, 0, length_includes_head=True,
            head_width=np.abs(y/50), head_length=x/10, color=aftcol) 
    ax.arrow(x, y, 0, -np.abs(y)/50, length_includes_head=True,
            head_width=np.abs(x/5), head_length=np.abs(y/100), color=aftcol) 
    #ax.text(x, y, 'PTF11agg', c=aftcol, va='top') 

    x = 0.2
    y = -25.9
    ax.scatter(x, y, c=aftcol, edgecolor='k', zorder=10) 
    ax.arrow(x, y, -x/3, 0, length_includes_head=True,
            head_width=np.abs(y/50), head_length=x/10, color=aftcol) 
    ax.arrow(x, y, 0, -np.abs(y)/50, length_includes_head=True,
            head_width=np.abs(x/5), head_length=np.abs(y/100), color=aftcol) 
    #ax.text(x, y, 'ZTF20aajnksq', c=aftcol, va='top') 

    x = 0.05
    y = -27.3
    ax.scatter(x, y, c=aftcol, edgecolor='k', zorder=10) 
    ax.arrow(x, y, -x/3, 0, length_includes_head=True,
            head_width=np.abs(y/50), head_length=x/10, color=aftcol) 
    ax.arrow(x, y, 0, -np.abs(y)/50, length_includes_head=True,
            head_width=np.abs(x/5), head_length=np.abs(y/100), color=aftcol) 
    #ax.text(x, y, 'ZTF21aaeyldq', c=aftcol, va='top') 

    x = 0.1
    y = -25
    ax.scatter(x, y, c=aftcol, edgecolor='k', zorder=10) 
    ax.arrow(x, y, -x/3, 0, length_includes_head=True,
            head_width=np.abs(y/50), head_length=x/10, color=aftcol) 
    ax.arrow(x, y, 0, -np.abs(y)/50, length_includes_head=True,
            head_width=np.abs(x/5), head_length=np.abs(y/100), color=aftcol) 
    #ax.text(x, y, 'ZTF21aayokph', c=aftcol, va='top') 

    ax.text(0.08, -27.2, 'Orphan Afterglows', color=aftcol, fontweight='bold')


if __name__=="__main__":
    fig,ax = plt.subplots(1,1,figsize=(5,4))#, sharey=True)

    # Plot BTS sources
    plot_bts(ax)

    # Plot ZTF sources
    plot_ztf(ax, background=False, shrink=2, text=True)

    # Plot CSS161010
    ax.errorbar(
            5.5, -21.5, 
            label=None, mfc=cowcol, mec=cowcol,
            c=cowcol, fmt='D', ms=5)

    # Plot AT2020mrf
    ax.errorbar(7.1, -20, 
            label=None, mfc=cowcol, mec=cowcol,
            c=cowcol, fmt='D', ms=5)

    # Add SNLS
    plot_snls(ax)

    # Plot AT2022tsd
    # 19.28 is the ext corr peak in g band
    # 19.74 is the value in r-band, which is g-band rest-frame
    Mpeak = 19.74-5*np.log10(vals.dL_mpc*1E6/10)+2.5*np.log10(1+vals.z)
    dur = 7.12
    edur = 2.57
    ax.errorbar(dur, Mpeak, xerr=edur, yerr=0.09, label=None, 
                c=cowcol, mec='k',
                fmt='D', ms=8, zorder=1000, lw=2)
    ax.text(
            dur, Mpeak*1.01, 'AT2022tsd', va='bottom', ha='left', 
            color=cowcol, fontweight='bold')

    ax.set_ylabel("Peak absolute mag.", fontsize=14)
    ax.set_ylim(-15,-23)
    
    # Luminosity axis
    ax2 = ax.twinx()
    # use the g-band effective wavelength: 4722.74 AA
    y_f = lambda y_i: 4E33*10**((y_i-4.77)/(-2.5))
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.plot([],[])
    ax2.set_yscale('log')
    ax2.tick_params(axis='both', labelsize=12)
    ax2.set_ylabel(
            r"Peak $\nu L_\nu$ (erg s$^{-1}$)", fontsize=14, rotation=270, 
            labelpad=15.0)

    ax.set_xlim(2.3,150)
    ax.set_xscale('log')
    ax.set_xlabel("Time above half-max (rest-frame days)", fontsize=14)
    ax.tick_params(axis='both', labelsize=12)

    # Plot the Mej = Mni line
    # Plot the Mni=Mej limit
    tpeak = np.linspace(1, 300)
    Mej = calc_Mej(tpeak)
    Lpeak = calc_Lpeak(Mej, tpeak)
    ax2.plot(tpeak, Lpeak, c='k', ls='--')
    ax.text(2.3, -16.8, r'$M_{\mathrm{Ni}}<M_{\mathrm{ej}}$', rotation=42)
    ax.text(2.3, -17.5, r'$M_{\mathrm{Ni}}>M_{\mathrm{ej}}$', rotation=42)

    plt.tight_layout()
    plt.show()
    #plt.savefig("lum_time_optical.png", dpi=200)
    #plt.close()
