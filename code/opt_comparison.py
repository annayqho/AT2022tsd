""" Compare the optical LC and spectrum to other classes of 
extragalactic transients.
Classes: 
    - CC SN
    - LFBOTs
    - TDEs
    - LGRBs 
    - LLGRBs 
"""


import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code)")
import vals
from get_opt import *


# Initialize a 2-panel figure
fig,axarr = plt.subplots(1,2,figsize=(8,3))

# Optical LC panel
ax = axarr[0]

# Plot the LC of the baseline transient
opt = get_full_opt()
# detections of the baseline transient
choose = np.logical_and(opt['isflare']==False, opt['mjdstart']<59856.4)
opt_transient = opt[choose]

# Plot
ms = ['*', 's', 'o', 'D']
cs = [vals.uc, vals.gc, vals.rc, vals.ic]
for i,b in enumerate(np.array(['u', 'g', 'r', 'i'])):
    choose = np.logical_and(opt_transient['flt']==b, opt_transient['emag']<99)
    dt = Time(opt_transient[choose]['mjdstart'].values, format='mjd').jd-vals.t0
    mag = opt_transient[choose]['mag'].values-vals.ext[b]
    emag = opt_transient[choose]['emag'].values
    ax.errorbar(dt, mag, emag, fmt=ms[i], color=cs[i])

# Add the extinction-corrected Keck data
t = 59871.44599
dt = Time(t,format='mjd').jd-vals.t0

g = 1.85
eg = 0.14
mg = -2.5*np.log10(g*1E-6)+8.90-vals.ext['g']
emg = (2.5/np.log(10)) * (eg/g)
ax.errorbar(dt,mg,emg,fmt='s',c=vals.gc)

i = 2.75
ei = 0.15
mi = -2.5*np.log10(i*1E-6)+8.90-vals.ext['i']
emi = (2.5/np.log(10)) * (ei/i)
ax.errorbar(dt,mi,emi,fmt='*',c=vals.ic)

ax.invert_yaxis()

# Optical spectrum panel
ax = axarr[1]

plt.show()




