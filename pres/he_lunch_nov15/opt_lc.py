""" Plot the optical light curve """
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.cosmology import Planck15


def plot_lc(ax,f,t,m,em,markers=False,lines=False,lw=1):
    cols = ['Crimson', 'Aquamarine', 'Goldenrod']
    symbols = ['o', 's', '>']
    for i,filt in enumerate(np.array(['r', 'g', 'i'])):
        choose = np.logical_and(em<99, f==filt)
        if markers:
            ax.errorbar(t[choose], m[choose], em[choose], 
                         fmt=symbols[i], c=cols[i], mec='k')
        if lines:
            ax.plot(t[choose], m[choose], c=cols[i], lw=lw)
        first_det_t = t[choose][0]
        first_det_m = m[choose][0]
        choose = np.logical_and(em==99, f==filt)
        if sum(choose)>0:
            if markers:
                ax.scatter(t[choose], m[choose], marker='v', c=cols[i], 
                           edgecolor='k')
            last_lim_t = t[choose][-1]
            last_lim_m = m[choose][-1]
            if lines:
                ax.plot([last_lim_t, first_det_t], [last_lim_m, first_det_m],
                         ls='-', c=cols[i], lw=lw)


def lc_AT2018cow(ax):
    """ Plot the LC of AT2018cow at z=0.256 

    Distances: factor of 21, meaning factor of 441 in flux,
    meaning difference in mag of 6.6
    """
    dat = pd.read_fwf("at2018cow_photometry_table.dat")
    t0 = 58285
    x = Time(dat['MJD'].values.astype('float'), format='mjd').value
    dat['Emag'][dat['Emag']=='-']=99
    for i in np.arange(len(dat['ABMag'])):
        if '>' in dat['ABMag'][i]:
            dat['ABMag'][i] = dat['ABMag'][i][1:]
    M = dat['ABMag'].values.astype(float)-Planck15.distmod(z=0.0141).value
    plot_lc(ax,
            dat['Filt'].values, (x-t0)*1.256, dat['ABMag'].values.astype(float)+6.6, 
            dat['Emag'].values.astype(float), lines=True, lw=1)
    ax.text(2, 23, '18cow at z=0.256', rotation=-35, fontsize=8)


def lc_AT2022tsd(ax):
    """ LC of our object """
    t0 = Time("2022-09-07")
    dat = pd.read_csv("opt.txt")

    x = Time(dat['Date'].values.astype(str), format='isot')-t0
    plot_lc(ax,
            dat['Filt'].values, x.value+2.5, 
            dat['Mag'].values,
            dat['Emag'].values, markers=True, lines=False, lw=2)
    ax.text(23, 22.5,'AT2022tsd', rotation=-25)


def lc_sn2002ap(ax):
    """ LC of SN2002ap
    z = 0.002108 """
    dratio = (Planck15.luminosity_distance(z=0.002108).value/Planck15.luminosity_distance(z=0.256).value)
    fratio = dratio**2
    moffset = -2.5*np.log10(fratio)

    dat = pd.read_fwf("LickLC02ap.dat", comment='#')
    t = dat['dJD'].values
    dt = (t-t[0])*1.256

    # Plot g-band
    apB = dat['B'].values.astype(float)
    apV = dat['V'].values.astype(float)
    g = ((apB-0.14) + (apV-0.02))/2.
    ax.plot(dt, g+moffset, ls=':', lw=0.7, c='Aquamarine')

    # Plot r-band
    apR = dat['Rc'].values
    apR = apR[0:-1].astype(float)
    R = apR+0.17
    ax.plot(dt[:-1], R+moffset, ls=':', lw=0.7, c='Crimson')

    # Plot i-band
    apI = dat['Ic'].values
    apI = apI[0:-1].astype(float)
    i = apI+0.43
    ax.plot(dt[:-1], i+moffset, ls=':', lw=0.7, c='Goldenrod')

    # Text
    ax.text(2, 25,'SN2002ap at z=0.256', rotation=25, fontsize=8)


if __name__=="__main__":
    fig,ax = plt.subplots(1,1,figsize=(6,4))
    lc_AT2022tsd(ax)
    lc_AT2018cow(ax)
    lc_sn2002ap(ax)

    # Make legend
    ax.text(0.05,0.95,'g',color='Aquamarine', transform=ax.transAxes)
    ax.text(0.07,0.95,'r',color='Crimson', transform=ax.transAxes)
    ax.text(0.09,0.95,'i',color='Goldenrod', transform=ax.transAxes)

    plt.gca().invert_yaxis()

    #plt.legend(loc='upper right')
    plt.xlabel("Days Since Explosion", fontsize=10)
    plt.ylabel("Absolute Mag", fontsize=10)
    plt.xlim(-4,43)
    plt.ylim(27,19)

    ax2 = ax.twinx()
    ax2.set_ylabel(
            r"Absolute Magnitude", fontsize=10, rotation=270,
            labelpad=15.0)
    # use the g-band effective wavelength: 4722.74 AA
    y_f = lambda y_i: y_i-Planck15.distmod(z=0.256).value
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.plot([],[])
    ax2.tick_params(axis='both', labelsize=10)

    #ax.axvline(x=30, ls='-', c='grey', lw=0.5)
    #ax.text(30, 21, 'HST DD Submitted', rotation=-90,fontsize=6)
    #ax.axvline(x=33, ls='-', c='grey', lw=0.5)
    #ax.text(33, 21.5, 'P200 Obs. (too faint)', rotation=-90,fontsize=6)
    #ax.axvline(x=33, ls='-', c='grey', lw=0.5)
    ax.axvspan(42,43,color='lightgrey')
    ax.text(41.5, 19.1, 'Keck g+I', rotation=-90,fontsize=8, va='top')
    #ax.axvspan(42,47,color='lightgrey')
    #ax.text(44, 19.1, 'VLT Epoch 2 (Oct 19-24)', rotation=-90,fontsize=6, va='top')
    #ax.axvspan(50,55,color='lightgrey')
    #ax.text(52, 19.1, 'VLT Epoch 3 (Oct 27-Nov 1)', rotation=-90,fontsize=6, va='top')
    #ax.axvline(x=43, ls='-', c='grey', lw=0.5)
    #ax.text(43, 21.5, 'VLT Epoch 2', rotation=-90,fontsize=6)
    #ax.axvspan(62,69,color='lightgrey')
    #ax.text(61, 21, 'VLT Epoch 3', fontsize=6)

    #plt.tight_layout()
    #plt.show()
    plt.savefig("opt_lc.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()

