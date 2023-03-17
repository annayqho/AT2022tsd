""" 
Based on the IPAC forced photometry light curve,
apply corrections, remove bad data points, and correct for MW extinction
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.time import Time


def get_lc():
    # Get the full IPAC forced phot light curve
    a = pd.read_table("ipac_forced_phot.txt", comment='#', delimiter=' ')

    # Get the full light curve
    jd = a['jd,'].values
    flux = a['forcediffimflux,'].values
    eflux = a['forcediffimfluxunc,'].values
    filt = a['filter,'].values
    maglim = a['diffmaglim,'].values
    zp = a['zpdiff,'].values
    programid = a['programid,'].values

    # Cobble together the LC
    SNT = 3
    SNU = 5 

    order = np.argsort(jd)
    mag = np.array([99]*len(flux)).astype(float)
    emag = np.array([99]*len(eflux)).astype(float)
    is_det = flux/eflux > SNT
    mag[is_det] = zp[is_det]-2.5*np.log10(flux[is_det])
    emag[is_det] = 1.0857*eflux[is_det]/flux[is_det]
    mag[~is_det] = zp[~is_det]-2.5*np.log10(SNU*eflux[~is_det])

    filt_final = np.copy(filt[order])
    filt_final[filt_final=='ZTF_g'] = 'g'
    filt_final[filt_final=='ZTF_r'] = 'r'
    filt_final[filt_final=='ZTF_i'] = 'i'

    # Final arrays for the light curve
    jd = jd[order]
    mag = mag[order]
    emag = emag[order]
    flux = flux[order]
    eflux = eflux[order]
    filt = filt_final

    return jd,filt,flux,eflux,mag,emag


def plot():
    fig,ax = plt.subplots(1,2,sharex=True)
    gdet = np.logical_and(filt=='g', emag<99)
    rdet = np.logical_and(filt=='r', emag<99)
    idet = np.logical_and(filt=='i', emag<99)

    glim = np.logical_and(filt=='g', emag==99)
    rlim = np.logical_and(filt=='r', emag==99)
    ilim = np.logical_and(filt=='i', emag==99)

    t0 = Time.now().jd
    #jd[gdet][np.argmin(mag[gdet])]
    ax[0].errorbar(jd[gdet]-t0, mag[gdet], emag[gdet], fmt='o', c='Aquamarine')
    ax[0].errorbar(jd[rdet]-t0, mag[rdet], emag[rdet], fmt='o', c='Crimson')
    ax[0].errorbar(jd[idet]-t0, mag[idet], emag[idet], fmt='o', c='Goldenrod')
    ax[0].scatter(jd[glim]-t0, mag[glim], marker='v', c='Aquamarine')
    ax[0].scatter(jd[rlim]-t0, mag[rlim], marker='v', c='Crimson')
    ax[0].scatter(jd[ilim]-t0, mag[ilim], marker='v', c='Goldenrod')

    ax[0].invert_yaxis()
    ax[0].set_xlabel("Days Ago")
    ax[0].set_ylabel("Mag")

    g = filt=='g'
    r = filt=='r'
    i = filt=='i'
    ax[1].errorbar(jd[g]-t0, fujy[g], efujy[g], fmt='o', c='Aquamarine')
    ax[1].plot(jd[g]-t0, fujy[g], c='Aquamarine')
    ax[1].errorbar(jd[r]-t0, fujy[r], efujy[r], fmt='o', c='Crimson')
    ax[1].errorbar(jd[i]-t0, fujy[i], efujy[i], fmt='o', c='Goldenrod')
    ax[1].set_xlabel("Days Ago")
    ax[1].set_ylabel("Flux (uJy)")
    ax[1].set_xlim(-30,0)

    plt.tight_layout()
    plt.show()
