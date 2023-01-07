""" 
Based on the IPAC forced photometry light curve,
apply corrections, remove bad data points, and correct for MW extinction
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.time import Time

# Get the full IPAC forced phot light curve
a = pd.read_table("lc.txt", comment='#', delimiter=' ')

# Extract all field/ccd/filter combinations
fid = np.zeros(len(a['filter,']))
fid[a['filter,']=='ZTF_g'] = 1
fid[a['filter,']=='ZTF_r'] = 2
fid[a['filter,']=='ZTF_i'] = 3
fcqf = (a['field,']*10000+a['ccdid,']*100+a['qid,']*10+fid).astype(int).values
ufcqf = np.unique(fcqf)

# Bad pixel contamination
# leaving this out for now...
if a['procstatus'].dtype=='int64':
    bad_pix = a['procstatus'].values == 56
else:
    bad_pix = np.logical_or(a['procstatus'].values == '56',
        a['procstatus'].values == '56,57')

# Bad observing conditions
bad_conditions = np.logical_or.reduce((
        a['scisigpix,'].values > 5*np.median(a['scisigpix,'].values),
        a['zpmaginpsci,'].values > 5*np.median(a['zpmaginpsci,'].values),
        a['zpmaginpscirms,'].values > 5*np.median(a['zpmaginpscirms,'].values)))

# Get the full light curve
jd = a['jd,'].values
flux = a['forcediffimflux,'].values
eflux = a['forcediffimfluxunc,'].values
filt = a['filter,'].values
maglim = a['diffmaglim,'].values
zp = a['zpdiff,'].values

# some flux values are NaN
nan = np.isnan(flux)

# some chisq values are bad
bad_chisq = np.isnan(a['forcediffimchisq,'])

# Discard bad values
#good = np.logical_and.reduce((~bad_pix, ~bad_conditions, ~nan, ~bad_chisq))
good = np.logical_and.reduce((~bad_conditions, ~nan, ~bad_chisq))
jd = jd[good]
flux = flux[good]
eflux = eflux[good]
filt = filt[good]
maglim = maglim[good]
zp = zp[good]
fcqf = fcqf[good]
refjdend = a['refjdend,'].values[good]
ufcqf = np.unique(fcqf)

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

# Convert the flux and eflux into physical units, uJy
f0 = 10**(0.4*zp)
fratio = flux/f0
fujy = fratio * 3631 * 1E6
efujy = eflux * (fujy/flux)

filt_final = np.copy(filt[order])
filt_final[filt_final=='ZTF_g'] = 'g'
filt_final[filt_final=='ZTF_r'] = 'r'
filt_final[filt_final=='ZTF_i'] = 'i'

# Final arrays for the light curve
jd = jd[order]
fujy = fujy[order]
efujy = efujy[order]
mag = mag[order]
emag = emag[order]
filt = filt_final

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
ax[0].set_xlabel("Days Since Peak")
ax[0].set_ylabel("Mag")

g = filt=='g'
r = filt=='r'
i = filt=='i'
ax[1].errorbar(jd[g]-t0, fujy[g], efujy[g], fmt='o', c='Aquamarine')
ax[1].errorbar(jd[r]-t0, fujy[r], efujy[r], fmt='o', c='Crimson')
ax[1].errorbar(jd[i]-t0, fujy[i], efujy[i], fmt='o', c='Goldenrod')
ax[1].set_xlabel("Days Ago")
ax[1].set_ylabel("Flux (uJy)")
ax[1].set_xlim(-30,0)

plt.tight_layout()
plt.show()
