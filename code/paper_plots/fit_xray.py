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
            'cmr.freeze', 2, cmap_range=(0.2, 0.8), return_fmt='hex')[::-1]

    df = load_swift()

    t = Time(df['!MJD    '].values, format='mjd')
    #dt = (t.jd-vals.t0)/(1+vals.z)
    dt = (t.jd-vals.t0)/(1+vals.z)
    et = (df['T_+ve   '].values)/(1+vals.z)

    L = df['L'].values/1E43
    eL = df['Lpos'].values/1E43

    isdet = eL>0

    ax.errorbar(
            dt[isdet],L[isdet],xerr=et[isdet],yerr=eL[isdet],
            fmt='o',c=cols[1], lw=0.5)
    ax.scatter(
            dt[~isdet],L[~isdet],marker='o',
            edgecolor=cols[1],facecolor='white', zorder=10)
    for i in np.arange(len(dt[~isdet])):
        ax.arrow(dt[~isdet][i], L[~isdet][i], dx=0, dy=-1,
                 length_includes_head=True, head_width=dt[~isdet][i]/20,
                 head_length=0.5, color=cols[1], zorder=0)

    popt,pcov = curve_fit(func, dt[isdet], L[isdet], p0=[10,20,-2], 
                          sigma=eL[isdet], absolute_sigma=True)

    xvals = np.linspace(22,85)
    yvals = popt[0]*(xvals/popt[1])**(popt[2])
    ax.plot(xvals,yvals,label='$L_X\propto t^{-1.8}$', c=cols[0])

    print(popt[2],np.sqrt(pcov[2,2]))

    ax.legend()
    ax.set_xscale('log')
    ax.set_xlim(22,85)
    #plt.yscale('log')
    ax.set_xticks([22,30,40,60,80])
    ax.set_xticklabels([22,30,40,60,80])
    ax.set_xlabel("$\Delta t$ (rest-frame days)")
    ax.set_ylabel("$L_X$ ($10^{43}$ erg s$^{-1}$)")

    plt.tight_layout()
    #plt.show()
    plt.savefig("xray_fit.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()
