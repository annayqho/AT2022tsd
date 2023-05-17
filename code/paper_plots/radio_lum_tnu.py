""" Peak luminosity vs. peak frequency phase space for radio transients """

import sys
sys.path.append("..")
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
import vals

cow_col = vals.cow_col

smallsize=8
medsize=10
bigsize=12

squaresize = 50


def ujy_to_flux(ujy, z):
    d = Planck15.luminosity_distance(z=z).cgs.value
    return ujy*1E-6*1E-23*4*np.pi*d**2


def vel_lines(ax, x, v):
    """ Equation 16 from Chevalier 1998 
    
    Parameters
    ----------
    ax: axis to plot the stuff on
    x: the value of (dt / 1 d) * (nu_p / 5 GHz)
    v: velocity in units of c
    """
    xvals = np.linspace(1,3000)
    logy = (26) + \
            (19/9) * np.log10(v) + \
            (19/9) * np.log10(xvals)
    yvals = 10**logy
    ax.plot(xvals, yvals, ls='--', c='k', lw=0.5)
    rotangle = 65
    ax.text(
            x, 2E28, "$R/\Delta t = %sc$" %v, 
            fontsize=smallsize, rotation=rotangle,
            horizontalalignment='center', verticalalignment='top', c='grey')
    return yvals


def mdot_curves(ax, x, y, mdotv, label=True):
    """ 
    x: x-coordinate for where to put the text
    y: y-coordinate for where to put the text
    mdotv: Mdot divided by v, in units of (10^{-4} Msol/yr) / (1000 km/s)
    x: (dt/1 day) * (nu_p / 5 GHz)
    """
    xvals = np.linspace(1,3000)
    eps_B = 1/3
    logy = (19/4) * np.log10(0.00005/eps_B) - (19/4)*np.log10(mdotv) + \
            (2*19/4)*np.log10(xvals) 
    yvals = 1E26 * 10**logy
    ax.plot(xvals, yvals, ls=':', c='k', lw=0.5)
    rotangle = 84
    if label:
        ax.text(
                x, y, 
                "$\dot{M}/v = 10^{%s}$" %int(np.log10(mdotv)), 
                fontsize=smallsize, rotation=rotangle,
                horizontalalignment='left', verticalalignment='top', c='grey')
    return yvals


def density_curves(ax, x, ne):
    """ 
    Parameters
    ----------
    ax: axis to plot the stuff on
    x: the value of (dt / 1 d) * (nu_p / 5 GHz)
    v: density in units of parts per cm cubed
    """
    xvals = np.linspace(1,3000)
    logy = (19/22) * np.log10(79) + 26 - (19/22)*np.log10(ne) + \
            (19/22) * 2 * np.log10(xvals**2 / 22) # divide out 22 days
    yvals = 10**logy
    ax.plot(xvals, yvals, ls=':', c='k', lw=0.5)
    rotangle = 75 
    ax.text(
            x, 5E29, "$n_e = 10^{%s} \mathrm{cm}^{-3}$" %int(np.log10(ne)), 
            fontsize=smallsize, rotation=rotangle,
            horizontalalignment='left', verticalalignment='top')
    return yvals


def typeii(ax):
    # 88Z, 79C
    tnu = np.array([1253*5/5, 1400*1.4/5])
    lpeak = np.array([2.2E28, 4.3E27])
    names = ['88Z', '79C']

    for i,tnuval in enumerate(tnu):
        if i==0:
            label='SN II'
        else:
            label=None
        ax.scatter(
                tnuval, lpeak[i], marker='o', edgecolor=vals.sn_col, 
                facecolor='white', s=100, label=label)
        ax.text(
                tnuval*1.1, lpeak[i]*1.1, names[i], fontsize=smallsize,
                verticalalignment='bottom',
                horizontalalignment='left', color=vals.sn_col)


def first(ax):
    # FIRST transient
    tnu = 26*365*(0.3/5)
    dcm = Planck15.luminosity_distance(z=0.01957).cgs.value
    lpeak = 2.25*1E-3*1E-23*4*np.pi*dcm**2
    ax.scatter(tnu, lpeak, marker='X', c='k', s=50, label="RT")
    ax.text(tnu, lpeak*1.2, 'FIRST J1419', fontsize=smallsize,
            verticalalignment='bottom',
            horizontalalignment='center')

    # Dillon's transient
    tnu = 4*365*(5/5)
    dcm = Planck15.luminosity_distance(z=0.03470).cgs.value
    lpeak = 7*1E-3*1E-23*4*np.pi*dcm**2
    ax.scatter(tnu, lpeak, marker='X', c='k', s=50)
    ax.text(tnu/1.3, lpeak*1.1, 'VT 1210+4956', fontsize=smallsize,
            verticalalignment='bottom',
            horizontalalignment='center')


def ibc(ax):
    # 2003L
    tnu = (30)*(22.5/5)
    lpeak = 3.3E28
    ax.scatter(
            tnu, lpeak, marker='+', c=vals.sn_col, s=100,
            label="SNe Ibc")
    ax.text(
            tnu, lpeak/1.2, "2003L", fontsize=smallsize,
            verticalalignment='top',
            horizontalalignment='center', color=vals.sn_col)

    # 11qcj
    tnu = (100)*(5/5)
    lpeak = 7E28
    ax.scatter(
            tnu, lpeak, marker='+', c=vals.sn_col, s=100,
            label=None)
    ax.text(
            tnu/1.2, lpeak, "11qcj", fontsize=smallsize,
            verticalalignment='center',
            horizontalalignment='right', color=vals.sn_col)

    # 2007bg
    tnu = (55.9)*(8.46/5)
    lpeak = 4.1E28
    ax.scatter(
            tnu, lpeak, marker='+', c=vals.sn_col, s=100, label=None)

    # SN 2003bg
    tnu = (35)*(22.5/5)
    lpeak = 3.9E28
    ax.scatter(
            tnu, lpeak, marker='+', c=vals.sn_col, s=100, label=None)

    # SN 2009bb
    tnu = (20)*(6/5)
    lpeak = 3.6E28
    ax.scatter(
            tnu, lpeak, marker='+', c=vals.sn_col, s=100)


def lfbot(ax):
    col = cow_col
    m = 'D'
    s = 30

    # Koala
    dcm = Planck15.luminosity_distance(z=0.2714).cgs.value
    tnu = np.array([(81/1.2714)*(10/5), (343)*(1.5/5)])/1.2714
    nu = np.array([10, 5])*1E9
    lpeak = np.array([0.364, 0.089])*1E-3*1E-23*4*np.pi*dcm**2
    ax.scatter(tnu, lpeak, marker=m, c=col, s=s)
    ax.plot(tnu, lpeak, color=col, ls='-')
    ax.text(
            tnu[0]/1.2, lpeak[0], "$\Delta t$=64d", fontsize=smallsize,
            verticalalignment='center',
            horizontalalignment='right', c=col)
    ax.text(
            tnu[1]*1.2, lpeak[1], "$\Delta t$=343d", fontsize=smallsize,
            verticalalignment='center',
            horizontalalignment='left', c=col)
    ax.text(tnu[0], lpeak[0]*1.2, "AT2018lug", fontsize=smallsize,
            horizontalalignment='right', c=col)

    # CSS 161010
    tnu = np.array([69*(5.6/5), 357*0.63/5])/1.033
    nu = np.array([5.6, 0.63])*1E9
    dcm = Planck15.luminosity_distance(z=0.033).cgs.value
    lpeak = np.array([8.8E-3, 1.2E-3])*1E-23*4*np.pi*dcm**2
    ax.scatter(
            tnu, lpeak, marker=m, c=col, s=s, label="_none")
    ax.text(tnu[0], lpeak[0]*1.2, "CSS161010", fontsize=smallsize,
            horizontalalignment='right', color=col,
            verticalalignment='bottom')
    ax.plot(tnu, lpeak, color=col, ls='-')
    ax.text(
            tnu[0]/1.2, lpeak[0], "$\Delta t$=69d", fontsize=smallsize,
            verticalalignment='center',
            horizontalalignment='right', c=col)
    ax.text(
            tnu[-1], lpeak[-1]/1.2, "$\Delta t$=357d", fontsize=smallsize,
            verticalalignment='top',
            horizontalalignment='center', c=col)

    # AT2020xnd
    x1 = 58*21.6/5
    dcm = Planck15.luminosity_distance(z=0.2442).cgs.value
    y1 = (0.68*1E-3*1E-23*4*np.pi*dcm**2)
    ax.scatter(
            x1, y1, marker=m, s=s, facecolors=col, edgecolors=col)
    ax.text(
            x1, y1/1.2, "$\Delta t$=22d", fontsize=smallsize, 
            verticalalignment='top',
            horizontalalignment='center', c=col)
    ax.text(
            x1/1.1, y1*1.1, "AT2020xnd", fontsize=smallsize, 
            verticalalignment='center',
            horizontalalignment='right', color=col)

    # AT2022tsd
    x = [26*(250/5)]
    dcm = Planck15.luminosity_distance(z=0.2567).cgs.value
    yf = np.array([0.6])
    y = yf*1E-3*1E-23*4*np.pi*dcm**2
    ax.scatter(x, y, marker=m, c=col, s=s*2, edgecolors='k',
               facecolors=col)
    ax.text(
            x[0], y[0]/1.2, "$\Delta t$=26d", fontsize=smallsize, 
            verticalalignment='top',
            horizontalalignment='center', c=col)
    ax.text(
            x[0]/1.08, y[0]*1.2,"AT2022tsd",fontsize=medsize, 
            fontweight='bold',verticalalignment='bottom',
            horizontalalignment='center', color=col)

    # AT2018cow
    x1 = 22*100/5
    y1 = 4.4E29
    ax.scatter(
            x1, y1, marker=m, s=s, 
            facecolors=col, edgecolors=col)
    ax.text(
            22*100/7*1.2, 5.5E29, "AT2018cow", fontsize=smallsize, 
            verticalalignment='bottom',
            horizontalalignment='left', color=col)
    ax.text(
            x1, y1/1.2, "$\Delta t$=22d", fontsize=smallsize, 
            verticalalignment='top',
            horizontalalignment='left', c=col)
    x2 = 91*10/5
    y2 = 4.3E28
    ax.scatter(x2, y2, marker=m, s=s, facecolors=col, edgecolors=col,
               label='LFBOT')
    ax.text(
            x2*1.1, y2*1, "$\Delta t$=91d", fontsize=smallsize,
            verticalalignment='bottom',
            horizontalalignment='left', c=col)
    plt.arrow(x1,y1,x2-x1,y2-y1, color=col)


def tde(ax):
    # ASASSN14li
    tnu = (143)*(8.20/5)
    lpeak = 1.8E28
    ax.scatter(
            tnu, lpeak, marker='o', edgecolor=vals.tde_col, 
            facecolor=vals.tde_col, s=100,
            label='TDE')
    ax.text(
            tnu*1.2, lpeak/1.3, "ASASSN14li", fontsize=smallsize,
            verticalalignment='top',
            horizontalalignment='center', color=vals.tde_col)


def llgrb(ax):
    # SN 1998bw
    tnu = (10)*(10/5)
    lpeak = 8.2E28
    ax.scatter(
            tnu, lpeak, marker='s', edgecolor=vals.llgrb_col, s=squaresize,
            facecolor=vals.llgrb_col, label="LLGRB-SN")
    ax.text(
            tnu, lpeak*1.2, "1998bw", fontsize=smallsize,
            verticalalignment='bottom',
            horizontalalignment='center', color=vals.llgrb_col)

    # GRB 171205A
    tnu = (4.3)*(6/5)
    dgrb = Planck15.luminosity_distance(z=0.0368).cgs.value
    # 3 mJy at 6 GHz with the VLA; Laskar et al. 2017
    lpeak = 3E-3 * 1E-23 * 4 * np.pi * dgrb**2
    ax.scatter(
            tnu, lpeak, marker='s', edgecolor=vals.llgrb_col, s=squaresize,
            facecolor=vals.llgrb_col, label=None)
    ax.text(
            tnu, lpeak*1.2, "2017iuk", fontsize=smallsize,
            verticalalignment='bottom',
            horizontalalignment='center', color=vals.llgrb_col)

    # SN 2006aj
    tnu = (5)*(4/5)
    lpeak = 8.3E27
    ax.scatter(
            tnu, lpeak, marker='s', edgecolor=vals.llgrb_col, s=squaresize,
            facecolor=vals.llgrb_col, label=None)
    ax.text(
            tnu, lpeak/1.3, "2006aj", fontsize=smallsize,
            verticalalignment='top',
            horizontalalignment='center', color=vals.llgrb_col)

    # SN 2010bh
    tnu = (30)*(5/5)
    lpeak = 1.2E28
    ax.scatter(
            tnu, lpeak, marker='s', edgecolor=vals.llgrb_col, s=squaresize,
            facecolor=vals.llgrb_col, label=None)
    ax.text(
            tnu/1.1, lpeak/1.3, "2010bh", fontsize=smallsize,
            verticalalignment='top',
            horizontalalignment='center', color=vals.llgrb_col)





if __name__=="__main__":
    # initialize figure
    fig,ax = plt.subplots(1,1, figsize=(5,4.4))

    # Plot each class
    typeii(ax)
    first(ax)
    ibc(ax)
    tde(ax)
    llgrb(ax)
    lfbot(ax)

    #lumtnu(ax)

    # Plot the background curves
    y = mdot_curves(ax, 500, 1.5E29, 10)
    y = mdot_curves(ax, 40, 2E28, 0.1)
    y = mdot_curves(ax, 4.5, 5E28, 0.001)
    y = vel_lines(ax, 12, 1)
    y = vel_lines(ax, 110, 0.1)
    y = vel_lines(ax, 1100, 0.01)

    # Add a legend
    ax.legend(bbox_to_anchor=(-0.2, 1.1), loc='upper left',
            ncol=6, fontsize=medsize, handletextpad=0.1,
            columnspacing=0.5, borderpad=0.3)

    # Formatting
    ax.set_ylabel(
        r"$L_{\mathrm{radio, peak}}$ (erg s$^{-1}$ Hz$^{-1}}$)",
        fontsize=medsize)
    ax.set_xlim(2, 3000)
    ax.set_ylim(3E27, 3E30)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', labelsize=medsize)
    ax.set_xlabel(
        "$(\Delta t/1\,\mathrm{day})(\\nu_p/5\,\mathrm{GHz})$",
        fontsize=medsize)

    # make a twin axis
    # ax2 = ax.twinx()
    # ax2.set_ylabel(
    #         r"$U/R$ (erg/cm) $\qquad \epsilon_e=\epsilon_B=1/3$", 
    #         fontsize=medsize, rotation=270, labelpad=15.0)
    # y_f = lambda y_i: 10**((14/19)*(np.log10(y_i)+14.65))
    # ymin, ymax = ax.get_ylim()
    # ax2.set_ylim((y_f(ymin), y_f(ymax)))
    # ax2.plot([],[])
    # ax2.set_yscale('log')
    # ax2.tick_params(axis='both', labelsize=medsize)
    # ax2.set_xlim(2,3000)

    plt.tight_layout()
    #plt.show()
    plt.savefig("lum_tnu.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()
