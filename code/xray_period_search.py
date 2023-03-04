import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import LombScargle
from get_xray import *

#x,y,xerr,yerr = load_chandra_flares('26641')
x,y,xerr,yerr = load_chandra_flares('26644')

keep = yerr > 0
x = x[keep]
y = y[keep]
yerr = yerr[keep]
ls = LombScargle(x, y, yerr)
frequency, power = ls.autopower()
p_minute = 1/frequency

window = np.logical_and(p_minute>70, p_minute<123)
pmax = p_minute[window][np.argmax(power[window])]
plt.axvline(x=pmax)
level = ls.false_alarm_level(0.025, method='bootstrap')
plt.axhline(y=level)

plt.plot(p_minute, power)
#plt.errorbar(x / pmax % 1, y, yerr, fmt='o')
plt.xlim(0, max(x))
plt.show()
