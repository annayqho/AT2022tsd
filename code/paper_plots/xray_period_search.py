import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import LombScargle
from get_xray import *


def all():
    x = []
    y = []
    ey = []
    for i,oid in enumerate(np.array(['26641', '26642', '26643', '26644'])):
        t,xi,yi,exi,eyi = load_chandra_flares(oid)
        [x.append(val) for val in t]
        [y.append(val) for val in yi]
        [ey.append(val) for val in eyi]
    x = np.array(x)
    dt = (x-x[0])*24*60 # in minutes
    y = np.array(y)
    ey = np.array(ey)
    keep = ey > 0
    dt = dt[keep]
    y = y[keep]
    ey = ey[keep]
    ls = LombScargle(dt, y, ey)
    frequency, power = ls.autopower(method='slow', maximum_frequency=1)
    period = (1/frequency) # in minutes
    plt.plot(period, power, 'k')
    level = ls.false_alarm_level(0.025, method='bootstrap')
    plt.axhline(y=level, c='lightgrey')


def per_obsid():
    c = ['k', 'lightgrey', 'Crimson', 'Aquamarine']
    for i,oid in enumerate(np.array(['26641', '26642', '26643', '26644'])):
        t,x,y,xerr,yerr = load_chandra_flares(oid)

        keep = yerr > 0
        x = x[keep]*24*60 # perform in minutes
        y = y[keep]
        yerr = yerr[keep]
        ls = LombScargle(x, y, yerr)
        frequency, power = ls.autopower(method='slow', maximum_frequency=1,
                                        samples_per_peak=5)
        p_minute = (1/frequency)
        #plt.plot(p_minute, power, color=c[i])
        plt.plot(frequency, power, color=c[i])
        level = ls.false_alarm_level(0.025, method='bootstrap')
        plt.axhline(y=level, color=c[i])

        #plt.errorbar(x / pmax % 1, y, yerr, fmt='o')

    plt.show()
