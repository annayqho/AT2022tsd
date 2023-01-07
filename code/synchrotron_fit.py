""" Fitting routines for a synchrotron spectrum """

import numpy as np
from scipy.optimize import curve_fit


def self_abs(x, b):
    """ Fitting function for self-absorbed part """
    y = 10**(2.5*np.log10(x)+b)
    return y


def fit_self_abs(freq, flux, flux_err=[]):
    """ Fit the self-absorbed part of the spectrum """
    if len(flux_err) > 0:
        popt, pcov = curve_fit(
                self_abs, freq, flux,
                sigma = flux_err, absolute_sigma = True)
    else:       
        popt, pcov = curve_fit(
                self_abs, freq, flux)
    return popt, pcov

