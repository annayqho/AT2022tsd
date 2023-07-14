""" Fit the radio spectral index on day 40, for use in our equipartition
analysis """

from astropy.time import Time
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from get_radio_at2022tsd import *


def powlaw(x,A,beta):
    return A*x**beta

dat = get_data()
dt = dat['dt']
flux = dat['Flux']
eflux = dat['eFlux']
choose = np.logical_and(dt>40.2, dt<44.2)
x = dat['Freq_Obs'][choose].values*(1+vals.z) # Rest-frame frequency
order = np.argsort(x)
x = x[order]
y = flux[choose].values[order] / (1+vals.z) # Correction in fnu space
ey = eflux[choose].values[order]

popt, pcov = curve_fit(powlaw,x,y,p0=[0.5,-1],sigma=ey,absolute_sigma=True)

plt.errorbar(x,y,ey,fmt='o',c='k')
xmod = np.linspace(100,300)
ymod = powlaw(xmod,*popt)
plt.plot(xmod,ymod,c='lightgrey')
print(popt[1])
print(np.sqrt(pcov[1,1]))
plt.tight_layout()
plt.show()
