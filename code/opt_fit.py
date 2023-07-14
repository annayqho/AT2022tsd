"""
Fit the optical transient LC with a power law 
"""


import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.optimize import curve_fit
import vals
from get_opt import *


def func(x,y0,x0,alpha):
    """ Fitting function for a power law """
    return y0*(x/x0)**(alpha)


# Initialize a figure
fig,ax = plt.subplots(1,1,figsize=(4,3))

# Plot the LC of the baseline transient
opt = get_full_opt()
# detections of the baseline transient
choose = np.logical_and(opt['isflare']==False, opt['mjdstart']<59856.4)
opt_transient = opt[choose]

# Create arrays of the transient LC
x = []
y = []
ey = []
filt = []

# Construct arrays
for i,b in enumerate(np.array(['g', 'r', 'i'])):
    choose = np.logical_and(opt_transient['flt']==b, opt_transient['emag']<99)
    dt = Time(opt_transient[choose]['mjdstart'].values, format='mjd').jd-vals.t0
    [x.append(val) for val in dt]
    mag = opt_transient[choose]['mag'].values-vals.ext[b]
    [y.append(val) for val in mag]
    emag = opt_transient[choose]['emag'].values
    [ey.append(val) for val in emag]
    [filt.append(val) for val in [b]*len(mag)]

# Add the extinction-corrected Keck data
t = 59871.44599
dt = Time(t,format='mjd').jd-vals.t0

g = 1.85
eg = 0.14
mg = -2.5*np.log10(g*1E-6)+8.90-vals.ext['g']
emg = (2.5/np.log(10)) * (eg/g)
x.append(dt)
y.append(mg)
ey.append(emg)
filt.append('g')

i = 2.75
ei = 0.15
mi = -2.5*np.log10(i*1E-6)+8.90-vals.ext['i']
emi = (2.5/np.log(10)) * (ei/i)
x.append(dt)
y.append(mi)
ey.append(emi)
filt.append('i')

filt = np.array(filt)
x = np.array(x)/(1+vals.z) # rest frame
y = np.array(y)
ey = np.array(ey)

# Now, convert back to flux
yf = 10**((y-8.90)/(-2.5)) / 1E-6
eyf = yf * (ey / (2.5/np.log(10)))

# Now plot
ms = ['s', 'o', 'D']
cs = [vals.gc, vals.rc, vals.ic]

for i,b in enumerate(np.array(['g', 'r', 'i'])):
    choose = filt==b
    ax.errorbar(x[choose], yf[choose], eyf[choose], fmt=ms[i],c=cs[i],
                label="$%s$" %b)
    order = np.argsort(x[choose])
    xfit = x[choose][order]
    yfit = yf[choose][order]
    eyfit = eyf[choose][order]
    if b in ['g', 'r']:
        popt, pcov = curve_fit(func, xfit, yfit, sigma=eyfit, 
                               absolute_sigma=True, p0=[20,10,-1])
        xplt = np.linspace(1,50)
        yplt = func(xplt, *popt)
        plt.plot(xplt, yplt, c=cs[i], label=r'$\alpha=%s\pm%s$' %(
            np.round(popt[2], 2), np.round(np.sqrt(pcov[2,2]),2)))

        yplt = func(xplt, *[6, 22.3, -1.6])
        if i==0:
            plt.plot(
                xplt, yplt, c='grey', ls='--', 
                lw=0.5, label=r'$\alpha=-1.6$')

ax.legend(fontsize=8, ncol=2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(3,100)
ax.set_xlabel("Rest-frame Days")
ax.set_ylabel(r"$f_\nu$ ($\mu$Jy)")

plt.savefig("opt_lc_fit.png", dpi=200, bbox_inches='tight', pad_inches=0.1)
plt.close()



