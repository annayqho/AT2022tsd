""" Optical rise time vs. peak luminosity
with a Nickel-56 curve """
import numpy as np
import matplotlib.pyplot as plt


def calc_Lpeak(Mni, tpeak):
    """
    Peak luminosity as a function of nickel mass and rise time

    Mni: nickel mass in units of Msun
    """
    tauni = 8.8
    tauco = 113.6
    Lpeak = 2E43 * (Mni) * (3.9*np.exp(-tpeak/tauni) + \
          0.678*(np.exp(-tpeak/tauco)-np.exp(-tpeak/tauni)))
    return Lpeak


def calc_Mej(tpeak):
    """
    The ejecta mass corresponding to rise time (Eq 2 of Kasen et al. 2017)
    We're assuming kappa = 0.1, V9=1

    Mej: in units of Msun
    tpeak: in units of days
    """
    Mej = (tpeak/14.5)**2
    return Mej


if __name__=="__main__":
    tpeak = np.linspace(1, 300)
    Mej = calc_Mej(tpeak)
    Lpeak = calc_Lpeak(Mej, tpeak)
    plt.plot(tpeak, Lpeak)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1,300)
    plt.ylim(1E41, 1E45)
    plt.show()
