""" Peak luminosity vs. peak frequency phase space for radio transients """

import sys
sys.path.append("..")
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15

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
            x, 4E27, "$R/\Delta t = %sc$" %v, 
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
    #mdotv = mdotv_scaled * 1000 / 1E-4
    xvals = np.linspace(1,3000)
    eps_B = 1/3
    logy = (19/4) * np.log10(0.0005/eps_B) - (19/4)*np.log10(mdotv) + \
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


def lumtnu(ax):
    # FIRST transient
    tnu = 26*365*(0.3/5)
    dcm = Planck15.luminosity_distance(z=0.01957).cgs.value
    lpeak = 2.25*1E-3*1E-23*4*np.pi*dcm**2
    ax.scatter(tnu, lpeak, marker='X', c='k', s=50)
    ax.text(tnu, lpeak*1.2, 'FIRST J1419', fontsize=smallsize,
            verticalalignment='bottom',
            horizontalalignment='center')

    # 88Z
    tnu = (1253)*(5/5)
    lpeak = 2.2E28
    ax.scatter(
            tnu, lpeak, marker='o', edgecolor='k', facecolor='white', s=100,
            label='SN II')
    ax.text(
            tnu, lpeak/1.2, "88Z", fontsize=smallsize,
            verticalalignment='top',
            horizontalalignment='right')

    # 79C
    tnu = (1400)*(1.4/5)
    lpeak = 4.3E27
    ax.scatter(
            tnu, lpeak, marker='o', edgecolor='k', facecolor='white', s=100,
            label=None)
    ax.text(
            tnu, lpeak/1.2, "79C", fontsize=smallsize,
            verticalalignment='top',
            horizontalalignment='right')

    # 2003L
    tnu = (30)*(22.5/5)
    lpeak = 3.3E28
    ax.scatter(
            tnu, lpeak, marker='+', c='k', s=100,
            label="SNe Ibc")
    ax.text(
            tnu, lpeak/1.2, "2003L", fontsize=smallsize,
            verticalalignment='top',
            horizontalalignment='center')

    # Koala
    dcm = Planck15.luminosity_distance(z=0.2714).cgs.value
    tnu = np.array([(81/1.2714)*(10/5), (343)*(1.5/5)])/1.2714
    nu = np.array([10, 5])*1E9
    lpeak = np.array([0.364, 0.089])*1E-3*1E-23*4*np.pi*dcm**2
    col = '#003f5c'
    ax.scatter(tnu, lpeak, marker='D', c=col, s=100)
    ax.plot(tnu, lpeak, color=col, ls='-')
    ax.text(
            tnu[0]/1.2, lpeak[0], "$\Delta t$=64d", fontsize=smallsize,
            verticalalignment='center',
            horizontalalignment='right', c=col)
    ax.text(
            tnu[1]*1.2, lpeak[1], "$\Delta t$=343d", fontsize=smallsize,
            verticalalignment='center',
            horizontalalignment='left', c=col)
    ax.text(tnu[0], lpeak[0]*1.2, "Koala", fontsize=smallsize,
            horizontalalignment='center', c=col)

    # CSS 161010
    col = '#ffa600'
    tnu = np.array([69*(5.6/5), 357*0.63/5])/1.033
    nu = np.array([5.6, 0.63])*1E9
    dcm = Planck15.luminosity_distance(z=0.033).cgs.value
    lpeak = np.array([8.8E-3, 1.2E-3])*1E-23*4*np.pi*dcm**2
    ax.scatter(
            tnu, lpeak, marker='h', c=col, s=100, 
            label="_none")
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

    # 11qcj
    tnu = (100)*(5/5)
    lpeak = 7E28
    ax.scatter(
            tnu, lpeak, marker='+', c='k', s=100,
            label=None)
    ax.text(
            tnu/1.2, lpeak, "11qcj", fontsize=smallsize,
            verticalalignment='center',
            horizontalalignment='right')


    # 2007bg
    tnu = (55.9)*(8.46/5)
    lpeak = 4.1E28
    ax.scatter(
            tnu, lpeak, marker='+', c='k', s=100, label=None)
    #ax.text(
    #        tnu, lpeak*1.2, "2007bg", fontsize=smallsize,
    #        verticalalignment='bottom',
    #        horizontalalignment='center')

    # SN 2003bg
    tnu = (35)*(22.5/5)
    lpeak = 3.9E28
    ax.scatter(
            tnu, lpeak, marker='+', c='k', s=100, label=None)
    #ax.text(
    #        tnu, lpeak/1.1, "2003bg", fontsize=smallsize,
    #        verticalalignment='top',
    #        horizontalalignment='left')

    # SN 1998bw
    tnu = (10)*(10/5)
    lpeak = 8.2E28
    ax.scatter(
            tnu, lpeak, marker='s', edgecolor='k', s=squaresize,
            facecolor='k', label="LLGRB-SN")
    ax.text(
            tnu, lpeak*1.2, "1998bw", fontsize=smallsize,
            verticalalignment='bottom',
            horizontalalignment='center')

    # GRB 171205A
    tnu = (4.3)*(6/5)
    dgrb = Planck15.luminosity_distance(z=0.0368).cgs.value
    # 3 mJy at 6 GHz with the VLA; Laskar et al. 2017
    lpeak = 3E-3 * 1E-23 * 4 * np.pi * dgrb**2
    ax.scatter(
            tnu, lpeak, marker='s', edgecolor='k', s=squaresize,
            facecolor='k', label=None)
    ax.text(
            tnu, lpeak*1.2, "2017iuk", fontsize=smallsize,
            verticalalignment='bottom',
            horizontalalignment='center')

    # ASASSN14li
    tnu = (143)*(8.20/5)
    lpeak = 1.8E28
    ax.scatter(
            tnu, lpeak, marker='o', edgecolor='k', facecolor='black', s=100,
            label='TDE')
    ax.text(
            tnu, lpeak/1.3, "ASASSN14li", fontsize=smallsize,
            verticalalignment='top',
            horizontalalignment='center')
    
    # J1644+57
    # tnu = (22)*(80/5)
    # lpeak = 1.1E32
    # ax.scatter(
    #         tnu, lpeak, marker='o', edgecolor='k', facecolor='black', s=100,
    #         label='TDE')
    # ax.text(
    #         tnu, lpeak/1.3, "J1644+57 (22d)", fontsize=smallsize,
    #         verticalalignment='top',
    #         horizontalalignment='center')

    # SN 2009bb
    tnu = (20)*(6/5)
    lpeak = 3.6E28
    ax.scatter(
            tnu, lpeak, marker='+', c='k', s=100)
    #ax.text(
    #        tnu/1.2, lpeak, "2009bb", fontsize=smallsize,
    #        verticalalignment='center',
    #        horizontalalignment='right')

    # SN 2006aj
    tnu = (5)*(4/5)
    lpeak = 8.3E27
    ax.scatter(
            tnu, lpeak, marker='s', edgecolor='k', s=squaresize,
            facecolor='k', label=None)
    ax.text(
            tnu, lpeak/1.3, "2006aj", fontsize=smallsize,
            verticalalignment='top',
            horizontalalignment='center')

    # SN 2010bh
    tnu = (30)*(5/5)
    lpeak = 1.2E28
    ax.scatter(
            tnu, lpeak, marker='s', edgecolor='k', s=squaresize,
            facecolor='k', label=None)
    ax.text(
            tnu, lpeak/1.3, "2010bh", fontsize=smallsize,
            verticalalignment='top',
            horizontalalignment='center')

    # Lines
    y = vel_lines(ax, 5.5, 1)
    y = vel_lines(ax, 55, 0.1)
    y = vel_lines(ax, 550, 0.01)

    # AT2020xnd
    x1 = 58*21.6/5
    dcm = Planck15.luminosity_distance(z=0.2442).cgs.value
    y1 = (0.68*1E-3*1E-23*4*np.pi*dcm**2)
    col = '#ef5675'
    ax.scatter(
            x1, y1, marker='*', s=300, 
            facecolors=col, edgecolors=col)
    ax.text(
            x1, y1/1.2, "$\Delta t$=22d", fontsize=smallsize, 
            verticalalignment='top',
            horizontalalignment='left', c=col)
    ax.text(
            x1*1.2, y1, "AT2020xnd", fontsize=smallsize, 
            verticalalignment='center',
            horizontalalignment='left', color=col)

    # AT2022tsd
    x = 26*(250/5)
    dcm = Planck15.luminosity_distance(z=0.2567).cgs.value
    y = 0.6*1E-3*1E-23*4*np.pi*dcm**2
    ax.scatter(
            x1, y1, marker='*', s=300, 
            facecolors=col, edgecolors=col)
    ax.text(
            x1, y1/1.2, "$\Delta t$=26d", fontsize=smallsize, 
            verticalalignment='top',
            horizontalalignment='left', c=col)
    ax.text(
            x1*1.2, y1, "AT2022tsd", fontsize=smallsize, 
            verticalalignment='center',
            horizontalalignment='left', color=col)

    # AT2018cow
    col = '#7a5195'
    x1 = 22*100/5
    y1 = 4.4E29
    ax.scatter(
            x1, y1, marker='o', s=200, 
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
    ax.scatter(
            x2, y2, marker='o', s=50, 
            facecolors=col, edgecolors=col)
    ax.text(
            x2*1.1, y2*1, "$\Delta t$=91d", fontsize=smallsize,
            verticalalignment='bottom',
            horizontalalignment='left', c=col)
    plt.arrow(x1,y1,x2-x1,y2-y1, color=col)


    ax.set_xlim(2, 3000)
    ax.set_ylim(9E26, 3E30)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', labelsize=medsize)
    ax.set_xlabel(
        "$(\Delta t/1\,\mathrm{day})(\\nu_p/5\,\mathrm{GHz})$",
        fontsize=medsize)


fig,ax = plt.subplots(1,1, figsize=(5,5))
lumtnu(ax)
y = mdot_curves(ax, 550, 2.5E29, 100, label=True)
y = mdot_curves(ax, 58, 4E29, 1, label=False)
y = mdot_curves(ax, 5.9, 6.4E29, 0.01)
#y = mdot_curves(ax, 1800, 1E-4)
ax.set_ylabel(
    r"$L_{\mathrm{radio, peak}}$ (erg s$^{-1}$ Hz$^{-1}}$)",
    fontsize=medsize)
#ax.get_yaxis().set_visible(False)
ax.legend(bbox_to_anchor=(-0.1, 1.1), loc='upper left',
        ncol=4, fontsize=medsize, 
        columnspacing=0.01, borderpad=0.3)#, columnspacing=0.1)
#y = mdot_curves(ax, 700, 1E1)

# make a twin axis
ax2 = ax.twinx()
ax2.set_ylabel(
        r"$U/R$ (erg/cm) $\qquad \epsilon_e=\epsilon_B=1/3$", 
        fontsize=medsize, rotation=270, labelpad=15.0)
y_f = lambda y_i: 10**((14/19)*(np.log10(y_i)+14.65))
ymin, ymax = ax.get_ylim()
ax2.set_ylim((y_f(ymin), y_f(ymax)))
ax2.plot([],[])
ax2.set_yscale('log')
ax2.tick_params(axis='both', labelsize=medsize)
ax2.set_xlim(2,3000)


plt.tight_layout()


     
plt.show()
#plt.savefig("lum_tnu.png", dpi=300)
#plt.close()
