#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 15:36:04 2022

@author: yuhanyao
"""
import numpy as np
import astropy.io.ascii as asci
from astropy.time import Time
import astropy.constants as const


from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70., Om0=0.3)

import matplotlib
import matplotlib.pyplot as plt
fs= 12
matplotlib.rcParams['font.size']=fs
ms = 6
matplotlib.rcParams['lines.markersize']=ms


def add_cow(ax, color = "k"):
    tbx = cow_xrt_lc()
    
    t0 = np.hstack([tbx["phase"].data, np.array([78.1, 211.8])])
    L0 = np.hstack([tbx["L_XRT"].data, np.array([1e+40, 7.1e+39])])
    eL0_left = np.hstack([tbx["L_XRT_unc_left"].data, np.array([2e+39, 2e+39])])
    eL0_right = np.hstack([tbx["L_XRT_unc_right"].data, np.array([1e+39, 2e+39])])

    ax.errorbar(t0, L0, yerr = [eL0_left, eL0_right], 
                fmt = "D-", color = color, 
                label = "AT2018cow", markersize=3, zorder = 3, linewidth = 1)
    
    ix = (t0>34.5)&(t0<37.5)
    subL0 = L0[ix]
    subL0_unc = eL0_left[ix]/2 + eL0_right[ix]/2
    subweight= 1 / subL0_unc**2
    subl0 = np.average(subL0, weights = subweight)
    subl0_unc = 1 / np.sqrt(np.sum(subweight))
    print ("AT2018cow at ~ 36 days luminosity: %.4f +- %.4f e+42 erg/s"%(subl0/1e+42, subl0_unc/1e+42))
    
    

def cow_xrt_lc():
    # T0 for this burst is Swift MET=551097181.6 s, = 2018 Jun 19 at 10:32:40.470 UT
    t0 = Time("2018-06-19T10:32:40.470").mjd
    #pdir = '/Users/yuhanyao/Dropbox/Mac/Documents/GitHub/CowAnalogs/data/'
    tb = asci.read("curve_nosys.qdp")
    names = ["Time", "T_+ve", "T_-ve", "Rate", "Ratepos", "Rateneg", "ObsID"]
    for i in range(len(names)):
        name = names[i]
        tb.rename_column("col%d"%(i+1), name)
    tb["mjd"] = t0 + tb["Time"]/3600/24
    tb["f_XRT"] =  tb["Rate"]*4.26e-11 # This is from Section 2.2.1 of Ho+2019
    tb["f_XRT_unc_right"] =  tb["Ratepos"]*4.26e-11
    tb["f_XRT_unc_left"] =  -tb["Rateneg"]*4.26e-11
    d_cm = 60 * 1e+6 * const.pc.cgs.value # 60 Mpc
    colnames = tb.colnames[-3:]
    for i in range(len(colnames)):
        colname_old = colnames[i]
        colname_new = "L"+colname_old[1:]
        tb[colname_new] = tb[colname_old] * 4 * np.pi * d_cm**2 
    # the four epochs
    #mjds1 = Time(["2018-06-23T17:31:09","2018-06-24T11:01:09"]).mjd
    #mjds2 = Time(["2018-07-02T14:00:12","2018-07-03T07:30:00"]).mjd
    #mjds3 = Time(["2018-07-14T06:20:09","2018-07-14T23:35:00"]).mjd
    #mjds4 = Time(["2018-07-22T13:56:09","2018-07-23T09:06:09"]).mjd
    #print_obsid(tb, mjds1, width = 0.5)
    #print_obsid(tb, mjds2, width = 1)
    #print_obsid(tb, mjds3, width = 2.5)
    #print_obsid(tb, mjds4, width = 3)
    
    # xmm 1st epoch
    #t1 = Time("2018-07-23T00:09:22").mjd
    # xmm 2nd epoch
    #t2 = Time("2018-09-06T14:49:20").mjd
    # xmm 3rd epoch
    #t3 = Time("2019-01-20T03:56:45").mjd
        
    tb["phase"] = tb["mjd"]-58285.44
    return tb



def add_2020mrf(ax, color = "tab:red"):
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
    ax.errorbar([36], [LX_20mrf], yerr = [[LX_20mrf_unc_left], [LX_20mrf_unc_right]],
                marker='s', color = color, 
                label = "AT2020mrf", markersize = 7, linewidth = lw,
                elinewidth = 1.5,
                markeredgewidth = 0.5, ecolor = "k")
    ax.errorbar(ts[3:], Ls[3:], yerr = [eLs_left[3:], eLs_right[3:]], marker='s', color = color, 
                markersize = 7, linestyle = "none", linewidth = lw,
                elinewidth = 1.5,
                markeredgewidth = 0.5, ecolor = "k")
    
    # eRASS3 upper limit
    tsrg_e3_mjd = Time(["2021-01-15T04:59:44", "2021-01-27T05:00:14"]).mjd
    srg3_phase_obs = (tsrg_e3_mjd - t0)/(1+z)
    tsrg_e3 = np.median(srg3_phase_obs)
    LXupp_20mrf_e3 = 7.24e-14 * 4 * np.pi * (D * const.pc.cgs.value)**2
    ax.plot(tsrg_e3, LXupp_20mrf_e3, marker = "v", markersize=7, color = color, 
            linestyle = "none", zorder = 2, markeredgewidth = 0.5)
    
    # eRASS4 upper limit
    tsrg_e4_mjd = Time(["2021-07-15T00:00:00", "2021-07-27T00:00:00"]).mjd
    srg4_phase_obs = (tsrg_e4_mjd - t0)/(1+z)
    tsrg_e4 = np.median(srg4_phase_obs)
    LXupp_20mrf_e4 = 8.26e-14 * 4 * np.pi * (D * const.pc.cgs.value)**2
    ax.plot(tsrg_e4, LXupp_20mrf_e4, marker = "v", markersize=7, color = color, 
            linestyle = "none", zorder = 2,markeredgewidth = 0.5)
    
    
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
    
    ax.errorbar(t1[:3], L1[:3], [eL1_left[:3], eL1_right[:3]], fmt = "o-", markersize=4.5, 
                color = color, zorder = 4, label = "AT2020xnd")
    ax.plot(t1[3:5], L1[3:5], marker = "v", markersize=5, color = color, linestyle = "none", zorder = 2)

    
def add_css(ax, color = "cyan"):
    t2 = np.array([99, 130, 291]) 
    # although I don't quite believe in the Coppejans+2020 analysis.... I will use the numbers
    f2 = np.array([1.33e-15, 1.94e-15, 1.31e-15])
    ef2 = np.array([0.76e-15, 0.97e-15, np.nan])
    distance_cm = cosmo.luminosity_distance(0.034).value*1e+6 * const.pc.cgs.value
    L2 = f2 * (4*np.pi*distance_cm**2)
    eL2 = ef2 * (4*np.pi*distance_cm**2)

    ax.errorbar(t2[:2], L2[:2], eL2[:2], color = color, 
                marker= "D", fmt="o--", markersize=4.5, zorder = 4, label = "CSS161010")
    ax.plot(t2[2:], L2[2:], marker = "v", markersize=5, color = color, linestyle = "none", zorder = 2)
       

def add_22tsd(ax, color = "red"):
    mjd1 = Time(
            ["2022-10-04T09:17:00", "2022-10-06T14:53:02", 
            "2022-10-08T02:15:00", '2022-10-09T05:05:00', 
             '2022-10-10T09:34:00', '2022-10-17T00:00:00',
             '2022-10-25T00:00:00', '2022-10-27T00:00:00', 
             '2022-11-04T00:00:00', '2022-11-06T00:00:00'], format='isot').mjd
    z = 0.256
    t0 = Time("2022-09-07").mjd
    ts = (mjd1-t0)/(1+z)
    Ls = np.array([1.02E44, 8.08E43, 8.13E43, 9.77E43, 5.49E43, 4.75E43, 8.57E43, 7E43, 2.4E43, 3*8.57E43/5])
    eLs = np.array([2.01E43, 1.72E43, 2.16E43, 2.98E43, 2.81E43, 1E42, 3.02E43, 0.1*7E43, 1E42, (0.0005/0.003)*3*8.57E43/5])
    ax.errorbar(ts, Ls, Ls*0.3, color = color, marker = "*", markersize = 15,
            markeredgecolor = "k", label = "AT2022tsd")


def custom_ax(ax):
    ax.semilogx()
    ax.semilogy()
    ax.set_xlabel("Rest-frame Days Since Explosion Epoch")
    ax.set_ylabel(r"$L_{\rm X}$"+" (erg/s)")
    
    z = 0.256
    distance_cm = cosmo.luminosity_distance(z).value*1e+6 * const.pc.cgs.value
    scale = 4*np.pi*distance_cm**2
    
    xstart = 22
    xmax = 350
    
    #xpos = 4
    #sens = 1.954E-12  * scale
    #color_nicer = "yellowgreen"
    #ax.plot([xstart, xmax], [sens, sens], linestyle = "-.", color = color_nicer)
    #ax.text(xpos, sens*0.8, "$NICER$ 1ks", color = color_nicer)
    
    xpos = 3.3
    sens = 9.5e-14* scale
    color_xrt = "yellowgreen"
    ax.plot([xstart, xmax], [sens, sens], linestyle = "-.", color = color_xrt)
    #ax.text(xpos, sens*0.8, "$Swift$ 5ks (9.5e-14 erg"+r"$\, \rm cm^{-2}\, s^{-1}$"+")", color = color_xrt)
    
    xpos = 2.3
    
    
    sens = 7.3e-15* scale
    color_cxo = "limegreen"
    #ax.plot([xstart, xmax], [sens, sens], linestyle = "-.", color = color_cxo)
    #ax.text(xpos, sens*0.8, "CXO 10ks (7.3e-15 erg"+r"$\, \rm cm^{-2}\, s^{-1}$"+")", color = color_cxo)
    
    """
    sens = 3.4e-15* scale
    color_cxo = "seagreen"
    ax.plot([xstart, xmax], [sens, sens], linestyle = "-.", color = color_cxo)
    ax.text(xpos, sens*1.05, "CXO 20ks (3.4e-15 cgs)", color = color_cxo)
    
    """
    sens = 2.2e-15* scale
    color_cxo = "darkgreen"
    #ax.plot([xstart, xmax], [sens, sens], linestyle = "-.", color = color_cxo)
    #ax.text(xpos, sens*0.8, "CXO 40ks (2.2e-15 erg"+r"$\, \rm cm^{-2}\, s^{-1}$"+")", color = color_cxo)
    
    
    xmin = 2.2
    xmax = 500
    ax.set_xlim(xmin, xmax)
    
    alphas = np.array([1.5, 2, 3])
    myrot = -24
    for i in range(len(alphas)):
        myalpha = alphas[i]
        xs = np.logspace(1.2, 2)
        if myalpha==1.5:
            #xs = np.logspace(1.2, 2.8)
            myrot = -15
            multi = 0.6
        elif myalpha == 2:
            #xs = np.logspace(1.2, 2.7)
            myrot = -25
            multi = 0.45
        elif myalpha == 3:
            #xs = np.logspace(1.2, 2.5)
            myrot = -35
            multi = 0.3
        t1 = 28
        L1 = 5E43
        A = L1 / (t1**(-1*myalpha))
        ys = A * xs**(-1*myalpha)
    
        ax.plot(xs, ys, color = "red", linestyle = ":", linewidth = 1.5)
    
        ax.text(xs[-1], ys[-1]*multi, r"$L \propto t^{-%.1f}$"%myalpha, 
                color = "r", fontsize = fs,
                zorder = 10, rotation =myrot, weight="bold")
    
    # upper x axis ticks
    axi1 = ax.twiny()
    axi1.tick_params(which = 'major', length = 4, direction = "in")
    axi1.tick_params(which = 'minor', length = 2, direction = "in")
    mjds = np.array([Time("2022-10-01").mjd,
                     Time("2022-11-01").mjd,
                     Time("2022-12-01").mjd,
                     Time("2023-01-01").mjd,
                     #Time("2023-02-01").mjd,
                     Time("2023-03-01").mjd,
                     #Time("2023-04-01").mjd,
                     Time("2023-06-01").mjd,
                     Time("2024-01-01").mjd])
    t0 = Time("2022-09-07").mjd
    x_months = (mjds - t0) / (1+z)
    
    axi1.set_xlim(xmin, xmax)
    axi1.semilogx()
    axi1.set_xticks(x_months)
    axi1.set_xticklabels(["2022-Oct",
                          "2022-Nov","2022-Dec",
                          "2023-Jan",#"2023-Feb",
                          "2023-Mar", #"2023-Apr", 
                          "2023-Jun", "2024-Jan"
                          ],
                         rotation = 80)
    axi1.set_xticks([], minor = True)
    axi1.set_xticklabels([], minor = True)
    


def make_figure1():
    plt.figure(figsize = (8.5, 5.5))
    ax = plt.subplot(111)
    add_cow(ax, color = "k")
    add_2020mrf(ax, color = "tab:blue")
    add_20xnd(ax, color = "tab:orange")
    add_css(ax, color = "tab:brown")
    add_22tsd(ax, color = "red")
    
    custom_ax(ax)
    ax.legend(loc = "lower left")

    #ax.axvline(x=60/(1.256), lw=0.5, c='red')#, label='Requested XRT Epochs')
    #ax.axvline(x=70/(1.256), lw=0.5, c='red')#, label='Requested XRT Epochs')
    #ax.axvline(x=42/(1.256), lw=0.5, c='red')#, label='Requested XRT Epochs')
    #ax.axvline(x=32, lw=0.5, c='red')#, label='Requested XRT Epochs')
    #ax.axvline(x=33, lw=0.5, c='red')
    #ax.axvline(x=36, lw=0.5, c='red')
    #ax.axvline(x=41, lw=0.5, c='red')

    plt.tight_layout(rect = (-0.02, -0.02, 1.02, 1.02))
    
    #plt.show()
    #plt.savefig("xlc_cxo.pdf")
    plt.savefig("xlc_xrt.png")
    plt.close()
    
    
if __name__ == '__main__':
    make_figure1()
