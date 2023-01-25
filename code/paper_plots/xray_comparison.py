""" Compare X-ray light curve to other classes of objects """

import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.table import Table
import astropy.constants as const
from astropy.cosmology import Planck18
import astropy.io.ascii as asci
import cmasher as cmr
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker

from fbot_xdata import cow_xrt_lc, add_SNeIbn_xlc, add_tde_lcs
from load_grb_xlc import add_grb_lcs, add_SNeIIn_xlc, add_xlc_sn1998bw, \
        add_xlc_sn2010dh, add_xlc_sn2006aj, add_xlc_sn2003dh, add_SLSNe_xlc


def add_at2022cmc(ax, color = "tab:darkgrey", ms=9):
    z = 1.1933
    a = pd.read_csv("at2022cmc_nicer.dat", delimiter=' ')
    ax.plot(a['x']*(1+z), a['L'], color=color, zorder=10)
    a = pd.read_csv("at2022cmc_xrt.dat", delimiter=' ')
    ax.plot(a['x']*(1+z), a['L'], color=color, zorder=10)


def add_cow(ax, color = "k"):
    tbx = cow_xrt_lc()
    t0 = np.hstack([tbx["phase"].data, np.array([78.1, 211.8])])
    L0 = np.hstack([tbx["L_XRT"].data, np.array([1e+40, 7.1e+39])])
    eL0_left = np.hstack([tbx["L_XRT_unc_left"].data, np.array([2e+39, 2e+39])])
    eL0_right = np.hstack([tbx["L_XRT_unc_right"].data, np.array([1e+39, 2e+39])])
    ax.scatter(t0, L0, marker='D', s=20, color=color, zorder=3)
    ax.plot(t0, L0, color=color, zorder=3)
    ax.scatter(0,0,marker='D', s=20, label='LFBOTs',color=color)


def add_css(ax, color = "cyan"):
    t2 = np.array([99, 130, 291])
    f2 = np.array([1.33e-15, 1.94e-15, 1.31e-15])
    ef2 = np.array([0.76e-15, 0.97e-15, np.nan])
    distance_cm = Planck18.luminosity_distance(0.034).value*1e+6 * const.pc.cgs.value
    L2 = f2 * (4*np.pi*distance_cm**2)
    eL2 = ef2 * (4*np.pi*distance_cm**2)
    ax.scatter(t2[:2], L2[:2], color = color, marker= "D", s=20)
    ax.plot(t2[:2], L2[:2], color = color)
    ax.scatter(t2[2:], L2[2:], marker = "v", color = color)
    ax.plot(t2[1:], L2[1:], ls='--', color=color)



def add_20xnd(ax, color = "orange"):
    tdis = 59132
    mjd1 = np.array([59157.8, 59163.8, 59179.1, 59207.2, 59316.6, 59317.1])
    distance_cm = Planck18.luminosity_distance(0.243).value*1e+6 * const.pc.cgs.value
    t1 = (mjd1-tdis) / (1+0.243)
    f1 = np.array([3.46e-14, 2.79e-14, 0.15e-14, 0.24e-14, 0.20e-14, 0.24e-14])
    ef1_right = np.array([0.96e-14, 0.75e-14, 0.17e-14, np.nan, np.nan, np.nan])
    ef1_left = np.array([1.27e-14, 0.67e-14, 0.11e-14, np.nan, np.nan, np.nan])
    L1 = f1 * (4*np.pi*distance_cm**2)
    eL1_right = ef1_right * (4*np.pi*distance_cm**2)
    eL1_left = ef1_left * (4*np.pi*distance_cm**2)

    ax.scatter(t1[:3], L1[:3], marker='D', color=color, zorder=4, s=20)
    ax.plot(t1[:3], L1[:3], color=color, zorder=4)


def add_2020mrf(ax, color = "tab:red", ms=9):
    z = 0.1353
    D = Planck18.luminosity_distance([z])[0].value * 1e+6 # in pc
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
    ax.scatter([36], [LX_20mrf], marker='D', color=color, s=20, zorder=4)
    ax.scatter(ts[3:], Ls[3:], marker='D', color = color, zorder=4, s=20)

    # eRASS3 upper limit
    tsrg_e3_mjd = Time(["2021-01-15T04:59:44", "2021-01-27T05:00:14"]).mjd
    srg3_phase_obs = (tsrg_e3_mjd - t0)#/(1+z)
    tsrg_e3 = np.median(srg3_phase_obs)
    LXupp_20mrf_e3 = 7.24e-14 * 4 * np.pi * (D * const.pc.cgs.value)**2
    ax.scatter(tsrg_e3, LXupp_20mrf_e3, marker = "v", color=color)

    # eRASS4 upper limit
    tsrg_e4_mjd = Time(["2021-07-15T00:00:00", "2021-07-27T00:00:00"]).mjd
    srg4_phase_obs = (tsrg_e4_mjd - t0)#/(1+z)
    tsrg_e4 = np.median(srg4_phase_obs)
    LXupp_20mrf_e4 = 8.26e-14 * 4 * np.pi * (D * const.pc.cgs.value)**2
    ax.scatter(tsrg_e4, LXupp_20mrf_e4, marker = "v", color=color)

    ax.plot(
            [36,tsrg_e3,ts[3],ts[4],tsrg_e4], 
            [LX_20mrf,LXupp_20mrf_e3,Ls[3],Ls[4],LXupp_20mrf_e4],
            ls='-',color=color)


def add_22tsd(ax, color = "red"):
    df = pd.read_csv("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/xray/AT2022tsd_XRT_binned_luminosity.qdp", delimiter=' ')
    mjd1 = df['!col1'].values/3600/24+59856.387276
    z = 0.256
    t0 = Time("2022-09-07").mjd
    ts = (mjd1-t0)#/(1+z)
    Ls = df['col4']
    uLs = df['col5']
    lLs = np.abs(df['col6'].values)
    ax.errorbar(ts, Ls, [lLs,uLs], color = color, marker = "D", markersize = 5,
            markeredgecolor = "k")



if __name__=="__main__":
    # Get colors. Same as mm comparison plot
    cols = cmr.take_cmap_colors('cmr.rainforest', 5, 
                                cmap_range=(0.1, 0.9), return_fmt='hex')[::-1]
    tde_col = cols[4]
    lgrb_col = cols[1]
    llgrb_col = cols[3]
    sn_col = cols[0]
    cow_col = cols[2]

    # Initialize figure 
    fig, ax = plt.subplots(1, 1, figsize=(3.5,6))

    # TDEs
    add_at2022cmc(ax, color=tde_col)
    add_tde_lcs(ax, color=tde_col)
    add_cow(ax, color = cow_col)
    add_css(ax, color = cow_col)
    add_20xnd(ax, color = cow_col)
    add_2020mrf(ax, color = cow_col)
    add_22tsd(ax, color= cow_col)

    #ax.scatter(30,1E44)

    #add_SNeIIn_xlc(ax)
    #add_SNeIbn_xlc(ax)
    #add_SLSNe_xlc(ax)

    add_grb_lcs(ax, color='lightgrey')

    # GRB-SN
    #add_xlc_sn1998bw(ax, color='lightgrey')
    #add_xlc_sn2010dh(ax, color='lightgrey')
    #add_xlc_sn2006aj(ax, color='lightgrey')
    #add_xlc_sn2003dh(ax, color='lightgrey')


    #add_afterglows(ax, color=cols[2])

    # ymax = 7e+47
    # ymin = 8e+38
    # xmin = 0.5
    # xmax = 600

    # # Region for CC SNe
    # xs = np.linspace(xmin, xmax)
    # yccmax = 5e+39
    # ax.fill_between(xs, ymin, yccmax, color = "gainsboro", zorder = 0, alpha = 0.8)
    # ax.arrow(20, yccmax, 0, -yccmax*0.25, color = "dimgray", head_length = yccmax*0.15,
    #          head_width = 2.0)
    # ax.text(15, yccmax*0.25, "Normal CCSNe", color = "dimgray")

    # ax.text(100, 2E46, 'J1644+57', c='k', fontsize=9)

    #ax.axvline(x=60,c='green',lw=0.5,ls='--')
    #ax.axvline(x=70,c='green',lw=0.5,ls='--')

    #ax.legend(loc = "lower left", ncol = 1, fontsize=9)
    #ax.text(20,1E44,'AT2022tsd',fontweight='bold',ha='right',color=cols[3])

    ax.set_xscale('log')
    ax.set_xlabel('$\Delta t_\mathrm{obs}$ (d)')
    ax.set_ylabel('$L_X$ (erg s$^{-1}$)')
    ax.set_yscale('log')

    plt.show()
    #plt.savefig('xray_lc.png', dpi=200, bbox_inches='tight', pad_inches=0.1)
    #plt.close()
