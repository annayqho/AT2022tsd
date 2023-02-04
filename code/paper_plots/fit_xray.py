""" Fit a power law to the X-ray light curve """

import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
import cmasher as cmr
import matplotlib.pyplot as plt
from astropy.time import Time
from scipy.optimize import curve_fit
from get_xray import *
import vals


def func(x,y0,x0,alpha):
    """ Fitting function for a power law """
    return y0*(x/x0)**(alpha)


if __name__=="__main__":
    fig,ax = plt.subplots(1,1,figsize=(4,3))

    cols = cmr.take_cmap_colors(
            'cmr.freeze', 3, cmap_range=(0.2, 0.8), return_fmt='hex')[::-1]

    # Load the Swift data
    df = load_swift()

    t = Time(df['!MJD    '].values, format='mjd')
    dt = (t.jd-vals.t0)/(1+vals.z)
    et = (df['T_+ve   '].values)/(1+vals.z)
    L = df['L'].values/1E43
    eL = df['Lpos'].values/1E43

    # Plot the Swift data
    c = cols[1]
    isdet = eL>0
    ax.errorbar(
            dt[isdet],L[isdet],xerr=et[isdet],yerr=eL[isdet],
            fmt='o',c=c, lw=0.5, label='Swift')
    ax.scatter(
            dt[~isdet],L[~isdet],marker='o',
            edgecolor=c,facecolor='white', zorder=10)
    for i in np.arange(len(dt[~isdet])):
        ax.arrow(dt[~isdet][i], L[~isdet][i], dx=0, dy=-L[~isdet][i]/4,
                 length_includes_head=True, head_width=dt[~isdet][i]/20,
                 head_length=L[~isdet][i]/8, color=c, zorder=0)

    # Get Chandra
    # WebPIMMS: the factor for going from 0.5-6 to 0.3-10 is 1.77
    factor = 1.77
    df = load_chandra()
    t = Time(df['MJD'].values, format='mjd')
    exp = (df['Exp'].values*1000)/86400 # days
    dt_start = (t.jd-vals.t0)/(1+vals.z)
    dt_dur = exp/(1+vals.z)
    dt_ch = dt_start + dt_dur/2 # halfway through
    e_dt_ch = dt_dur/2 # half the exposure
    L_ch = 1.77*df['L'].values/1E43
    uL_ch = 1.77*df['Lpos'].values/1E43
    lL_ch = 1.77*df['Lneg'].values/1E43

    # Plot Chandra
    c = cols[2]
    ax.errorbar(
            dt_ch,L_ch,xerr=e_dt_ch,yerr=[lL_ch,uL_ch], 
            fmt='s',c=c, lw=0.5, label='Chandra')

    # Concatenate
    dt = np.hstack((dt[isdet], dt_ch))
    L = np.hstack((L[isdet], L_ch))
    eL = np.hstack((eL[isdet], uL_ch))
    e_dt = np.hstack((et[isdet], e_dt_ch))

    # Fit 
    popt,pcov = curve_fit(func, dt, L, p0=[10,20,-2], 
                          sigma=eL, absolute_sigma=True)
    xvals = np.linspace(22,140)
    yvals = popt[0]*(xvals/popt[1])**(popt[2])
    ax.plot(xvals,yvals,label='$L_X\propto t^{-1.8}$', c=cols[0])
    print(popt[2],np.sqrt(pcov[2,2]))

    ax.legend()
    ax.set_xscale('log')
    ax.set_xlim(20,140)
    plt.yscale('log')
    ax.set_xticks([20,30,40,60,80, 110, 140])
    ax.set_xticklabels([20,30,40,60,80, 110, 140])
    ax.set_xlabel("$\Delta t$ (rest-frame days)")
    ax.set_ylabel("$L_X$ ($10^{43}$ erg s$^{-1}$)")

    plt.tight_layout()
    #plt.show()
    plt.savefig("xray_fit.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()
