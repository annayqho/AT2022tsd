#!/usr/bin/env python
# coding: utf-8

# # Figure 2: X-ray LC compilation (original code from Yuhan Yao) 

# In[1]:


import numpy as np
import pandas as pd

from astropy.time import Time
from astropy.table import Table
import astropy.constants as const
import astropy.io.ascii as asci

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70., Om0=0.3)
 
import sys
sys.path.append("..")
import vals

import cmasher as cmr

import matplotlib
import matplotlib.pyplot as plt
fs = 10
matplotlib.rcParams['font.size']=fs
from matplotlib import ticker
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'


# In[2]:


from fbot_xdata import cow_xrt_lc, add_SNeIbn_xlc, add_tde_lcs


# In[3]:


from load_grb_xlc import add_grb_lcs, add_SNeIIn_xlc, add_xlc_sn1998bw, add_xlc_sn2010dh, \
                        add_xlc_sn2006aj, add_xlc_sn2003dh, add_SLSNe_xlc #add_grb111209A_lc


# In[4]:


cols = cmr.take_cmap_colors(
        'cmr.rainforest', 5, cmap_range=(0.0, 0.85), return_fmt='hex')


# In[5]:


def add_cow(ax, color = "k"):
    tbx = cow_xrt_lc()
    
    t0 = np.hstack([tbx["phase"].data, np.array([78.1, 211.8])])
    L0 = np.hstack([tbx["L_XRT"].data, np.array([1e+40, 7.1e+39])])
    eL0_left = np.hstack([tbx["L_XRT_unc_left"].data, np.array([2e+39, 2e+39])])
    eL0_right = np.hstack([tbx["L_XRT_unc_right"].data, np.array([1e+39, 2e+39])])

    #ax.errorbar(t0, L0, yerr = [eL0_left, eL0_right], 
    #            fmt = "D-", color = color, markersize=3, zorder = 3, linewidth = 1)
    ax.scatter(t0, L0, marker='D', s=6, color=color, zorder=3)
    ax.plot(t0, L0, color=color, zorder=3, lw=0.5)
    ax.scatter(0,0,marker='D', s=20, label='LFBOTs',color=color)
    
#     ix = (t0>34.5)&(t0<37.5)
#     subL0 = L0[ix]
#     subL0_unc = eL0_left[ix]/2 + eL0_right[ix]/2
#     subweight= 1 / subL0_unc**2
#     subl0 = np.average(subL0, weights = subweight)
#     subl0_unc = 1 / np.sqrt(np.sum(subweight))
#     print ("AT2018cow at ~ 36 days luminosity: %.4f +- %.4f e+42 erg/s"%(subl0/1e+42, subl0_unc/1e+42))
    
    
def add_css(ax, color = "cyan"):
    t2 = np.array([99, 130, 291]) 
    # although I don't quite believe in the Coppejans+2020 analysis.... I will use the numbers
    f2 = np.array([1.33e-15, 1.94e-15, 1.31e-15])
    ef2 = np.array([0.76e-15, 0.97e-15, np.nan])
    distance_cm = cosmo.luminosity_distance(0.034).value*1e+6 * const.pc.cgs.value
    L2 = f2 * (4*np.pi*distance_cm**2)
    eL2 = ef2 * (4*np.pi*distance_cm**2)

    ax.scatter(t2[:2], L2[:2], color = color, marker= "D", s=10)
    ax.plot(t2[:2], L2[:2], color = color)
    ax.scatter(t2[2:], L2[2:], marker = "v", color = color, s=10)
    ax.plot(t2[1:], L2[1:], ls='--', color=color)



def add_20xnd(ax, color = "orange"):
    tdis = 59132
    mjd1 = np.array([59157.8, 59163.8, 59179.1, 59207.2, 59316.6, 59317.1])
    distance_cm = cosmo.luminosity_distance(0.243).value*1e+6 * const.pc.cgs.value
    t1 = (mjd1-tdis) / (1+0.243)
    f1 = np.array([3.46e-14, 2.79e-14, 0.15e-14, 0.24e-14, 0.20e-14, 0.24e-14])
    ef1_right = np.array([0.96e-14, 0.75e-14, 0.17e-14, np.nan, np.nan, np.nan])
    ef1_left = np.array([1.27e-14, 0.67e-14, 0.11e-14, np.nan, np.nan, np.nan])
    L1 = f1 * (4*np.pi*distance_cm**2)
    eL1_right = ef1_right * (4*np.pi*distance_cm**2)
    eL1_left = ef1_left * (4*np.pi*distance_cm**2)
    
    ax.scatter(t1[:3], L1[:3], marker='D', color=color, zorder=10, s=10, edgecolor='k')
    #ax.errorbar(t1[:3], L1[:3], [eL1_left[:3], eL1_right[:3]], fmt = "o-", markersize=4.5, 
    #            color = color, zorder = 4, label = "AT2020xnd")
    ax.plot(t1[:3], L1[:3], color=color, zorder=4)
    #ax.plot(t1[3:5], L1[3:5], marker = "v", markersize=5, color = color, linestyle = "none", zorder = 2)

    
    
    
def custom_ax(ax, ymax = 5e+43, ymin = 8e+38, 
              xmin = 2.5, add_phase = False, xmax = 350):
    ax.semilogx()
    ax.semilogy()

    ax.tick_params(which = 'major', length = 4, top=True, direction = "in", right = True, labelsize=10)
    ax.tick_params(which = 'minor', length = 2, top=True, direction = "in", right = True, labelsize=10)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    if add_phase:
        ax.fill_between(np.array([20, 90]), ymin, ymax, color = "silver", alpha = 0.5)
        ax.text(20*0.2, ymin*1.2, "Plateau Phase")
        ax.text(20*1.1, ymin*1.2, "Decline Phase")
        #ax.text(130*0.9, ymin*1.2, "Plateau Phase?")    
        
    lgymax = np.log10(ymax)
    lgymin = np.log10(ymin)
    
    majorys = []
    for yy in range(38, 50):
        if (yy > lgymin) & (yy < lgymax):
            YY = 10**yy
            YY = np.float(YY)
            majorys.append(YY)
    majorys = np.array(majorys)
    
    ax.set_yticks(majorys)
    
    minorys = []
    for yy in np.arange(38, 50, 0.1):
        #print (yy)
        if (yy > lgymin) & (yy < lgymax):
            YY = 10**yy
            YY = np.float(YY)
            minorys.append(YY)
    minorys = np.array(minorys)
        
    ax.set_yticks(minorys, minor=True)
    ax.yaxis.set_minor_formatter(ticker.NullFormatter()) 
    
    
def add_cxo_tick(ax, ymax):
    tt = np.array([6, 10, 20, 40, 80, 160])
    #tt_ = np.array([8, 16, 50, 100, 170, 260])
    ax.plot([tt, tt], [ymax, 1.2*ymax], color = "r", clip_on=False)
    #ax.plot([tt_, tt_], [ymax, 1.2*ymax], color = "darkorange", clip_on=False)
    for i in range(len(tt)):
        ax.text(tt[i]*0.95, 1.3*ymax, "$C$", clip_on=False, color = "r", fontsize = fs-2)
    #for i in range(len(tt_)):
    #    ax.text(tt_[i]*0.95, 1.3*ymax, "$C$", clip_on=False, color = "darkorange", fontsize = fs-2)


def add_2020mrf(ax, color = "tab:red", ms=9):
    z = 0.1353
    D = cosmo.luminosity_distance([z])[0].value * 1e+6 # in pc
    t0 = 59012  
    D_cm = D*const.pc.cgs.value

    #tsrg_e1_mjd = Time(["2020-01-19T10:41:10.661","2020-01-22T14:41:21.682"]).mjd
    #srg1_phase_obs = (tsrg_e1_mjd - t0)/(1+z)
    
    # eRASS2 detection
    tsrg_e2 = np.array([35, 36, 37])
    LX_20mrf = 3.9039e-13 * 4 * np.pi * D_cm**2
    LX_20mrf_unc_right = 4.786e-13 * 4 * np.pi * D_cm**2 - LX_20mrf
    LX_20mrf_unc_left = LX_20mrf - 3.110e-13 * 4 * np.pi * D_cm**2
    frac = np.array([0.05309, 0.32039, 0.05112])
    meancr = 0.10891
    lsrg_e2 = LX_20mrf / meancr * frac
    lsrg_e2_unc_right = LX_20mrf_unc_right/ meancr * frac
    lsrg_e2_unc_left = LX_20mrf_unc_left/ meancr * frac
    
    
    print ("AT2020mrf eRASS2 Luminosity = %.2f - %.2f + %.2f e+43 erg/s"%(LX_20mrf/1e+43, 
                                                                          LX_20mrf_unc_left/1e+43, 
                                                                          LX_20mrf_unc_right/1e+43))
    # Two Chandra obsIDs
    tc1 = 327.4
    tc2 = 328.2
    Lc1 = 4.00e-14 * 4 * np.pi * (D * const.pc.cgs.value)**2
    Lc2 = 1.57e-14 * 4 * np.pi * (D * const.pc.cgs.value)**2
    Lc1_unc_right = 0.68e-14 * 4 * np.pi * (D * const.pc.cgs.value)**2
    Lc1_unc_left = 1.24e-14 * 4 * np.pi * (D * const.pc.cgs.value)**2
    Lc2_unc_right = 0.27e-14 * 4 * np.pi * (D * const.pc.cgs.value)**2
    Lc2_unc_left = 0.49e-14 * 4 * np.pi * (D * const.pc.cgs.value)**2
    
    Lc_mean = np.mean([Lc1, Lc2])
    
    lnL_ratio = np.log(LX_20mrf / Lc_mean)
    lnt_ratio = np.log(328 / 36)
    lnL_ratio  / lnt_ratio

    ts = np.hstack([tsrg_e2, np.array([tc1, tc2]) ])
    Ls = np.hstack([lsrg_e2, np.array([Lc1, Lc2]) ])
    eLs_right = np.hstack([lsrg_e2_unc_right, np.array([Lc1_unc_right, Lc2_unc_right]) ])
    eLs_left = np.hstack([lsrg_e2_unc_left, np.array([Lc1_unc_left, Lc2_unc_left]) ])
    lw = 0.8
    #ax.errorbar(ts[:3], Ls[:3], yerr = [eLs_left[:3], eLs_right[:3]], marker='*', color = color, 
    #            label = "AT2020mrf", markersize = ms+2, linestyle = "-", linewidth = lw,
    #            elinewidth = 1.5,
    #            markeredgecolor = "k", markeredgewidth = 0.5, ecolor = "k")
    ax.scatter([36], [LX_20mrf], marker='D', color=color, s=10, zorder=4)
    ax.scatter(ts[3:], Ls[3:], marker='D', color = color, zorder=4, s=10)
    
    # eRASS3 upper limit
    tsrg_e3_mjd = Time(["2021-01-15T04:59:44", "2021-01-27T05:00:14"]).mjd
    srg3_phase_obs = (tsrg_e3_mjd - t0)#/(1+z)
    tsrg_e3 = np.median(srg3_phase_obs)
    LXupp_20mrf_e3 = 7.24e-14 * 4 * np.pi * (D * const.pc.cgs.value)**2
    #ax.scatter(tsrg_e3, LXupp_20mrf_e3, marker = "v", color=color, s=10)
    
    # eRASS4 upper limit
    tsrg_e4_mjd = Time(["2021-07-15T00:00:00", "2021-07-27T00:00:00"]).mjd
    srg4_phase_obs = (tsrg_e4_mjd - t0)#/(1+z)
    tsrg_e4 = np.median(srg4_phase_obs)
    LXupp_20mrf_e4 = 8.26e-14 * 4 * np.pi * (D * const.pc.cgs.value)**2
    #ax.scatter(tsrg_e4, LXupp_20mrf_e4, marker = "v", color=color, s=15)
    
    #ax.plot([36,tsrg_e3,ts[3],ts[4],tsrg_e4], 
    #        [LX_20mrf,LXupp_20mrf_e3,Ls[3],Ls[4],LXupp_20mrf_e4],
    #        ls='-',color=color, lw=0.5)
    ax.plot([36,ts[3],ts[4]], 
            [LX_20mrf,Ls[3],Ls[4]],
            ls='-',color=color, lw=0.5)
        
    ax.plot([ts[3],ts[4]], [Ls[3],Ls[4]],ls='-',color=color)


# In[6]:


def add_at2022cmc(ax, color = "tab:darkgrey", ms=9):
    z = 1.1933
    a = pd.read_csv("../../data/xray/at2022cmc_nicer.dat", delimiter=' ')
    ax.plot(a['x']*(1+z), a['L'], color=color, zorder=10)#, label='AT2022cmc')
    a = pd.read_csv("../../data/xray/at2022cmc_xrt.dat", delimiter=' ')
    ax.plot(a['x']*(1+z), a['L'], color=color, zorder=10, label='TDE')


# In[7]:


def add_afterglows(ax, color='k', ms=9):
    
    # ZTF19abvizsw
    z = 1.2596
    ax.plot(np.array([1.9,9.9]),[2.5E46, 3.8E45],c=color, zorder=2)
    ax.scatter(np.array([1.9,9.9]),[2.5E46, 3.8E45],marker='o', c=color, zorder=2)
    
    # ZTF20aajnksq
    z = 2.9
    ax.scatter(2.1,9.9E45,c=color,marker='o')
    ax.plot(np.array([2.1,3.4]), [9.9E45, 9.9E45], c=color, ls='--', zorder=2)
    ax.scatter(3.4,9.9E45,marker='v',c=color, zorder=2)
    
    # ZTF21aayokph
    z = 1.0624
    ax.scatter(1.94, 2.2E45, marker='o', c=color, zorder=2)
    ax.plot(np.array([1.94,4.39]), [2.2E45, 1.0E45], c=color, ls='--', zorder=2)
    ax.scatter(4.39, 1.0E45, marker='v', c=color, zorder=2)
    
    # ZTF21aaeyldq
    z = 2.5131
    ax.scatter(0.83, 1.7E46, marker='o', c=color, label="ZTF Afterglows", zorder=2)
    ax.plot(np.array([0.83, 3.80]), [1.7E46, 1.1E46], c=color, ls='--', zorder=2)
    ax.scatter(3.80,1.1E46,c=color,marker='v', zorder=2)


# In[8]:


from get_xray import load_swift, load_chandra, load_both
def add_22tsd(ax, color = "red"):
    dt, L, e_dt, eL = load_both()
    order = np.argsort(dt)
    ratio = L[order]/eL[order]
    choose = ratio > 1
    #ax.errorbar(dt[order][choose], L[order][choose], yerr=eL[order][choose], color = color, marker = "D", markersize = 5,
    #        markeredgecolor = "k", lw=0.5, zorder=200) 
    ax.errorbar(dt[order][choose], L[order][choose], yerr=eL[order][choose], color = color, marker = "D", markersize = 4,
            markeredgecolor = "k", lw=0.5, zorder=200)     


# In[9]:


def create_xray_panel(ax):
    """ Create a panel showing the X-ray LC comparison """
    
    cols = cmr.take_cmap_colors(
        'cmr.rainforest', 5, cmap_range=(0.1, 0.9), return_fmt='hex')[::-1]

    tde_col = vals.tde_col
    sn_col = vals.sn_col
    llgrb_col = vals.llgrb_col
    lgrb_col = vals.lgrb_col
    cow_col = vals.cow_col

    add_at2022cmc(ax, color=tde_col)

    add_cow(ax, color = cow_col)
    add_css(ax, color = cow_col) 
    add_20xnd(ax, color = cow_col)    
    add_2020mrf(ax, color = cow_col)
    add_22tsd(ax, color= cow_col)

    add_grb_lcs(ax, color=lgrb_col)

    # GRB-SN
    add_xlc_sn1998bw(ax, color=llgrb_col)
    add_xlc_sn2010dh(ax, color=llgrb_col)
    add_xlc_sn2006aj(ax, color=llgrb_col)
    add_xlc_sn2003dh(ax, color=llgrb_col)

    add_tde_lcs(ax, color=tde_col)

    ymax = 7e+47
    ymin = 8e+38
    xmin = 0.5
    xmax = 600

    xs = np.linspace(xmin, xmax)
    yccmax = 5e+39
    ax.fill_between(xs, ymin, yccmax, color = sn_col, zorder = 0, alpha = 0.8)
    ax.arrow(20, yccmax, 0, -yccmax*0.25, color = 'white', head_length = yccmax*0.15,
             head_width = 2.0)
    ax.text(15, yccmax*0.25, "Normal CCSNe", color = 'white')

    ax.legend(loc = "lower left", #bbox_to_anchor = (0, 0.2)
              ncol = 1, fontsize=8, handletextpad=0.1)

    ax.set_xlabel("$\Delta t_\mathrm{obs}$ (d)", fontsize=10)
    ax.set_ylabel("X-ray luminosity (erg s$^{-1}$)", fontsize=10)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', labelsize=10)
    ax.set_xlim(1,600)
    ax.set_ylim(1E39,1E48)

    ax.text(5, 1.5E43, 'AT2018cow', c=cow_col, ha='right', fontsize=8)
    ax.text(45, 0.8E41, 'AT2020xnd', c=cow_col, ha='right', fontsize=8)
    ax.text(600, 0.3E43, 'AT2020mrf', c=cow_col, ha='right', fontsize=8)
    ax.text(60, 0.6E44, 'AT2022tsd', c=cow_col, ha='left', fontsize=8, fontweight='bold')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




