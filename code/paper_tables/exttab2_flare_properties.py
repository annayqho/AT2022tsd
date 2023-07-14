""" Calculate the properties of each flare: duration, luminosity, etc """

import sys
sys.path.append("..")
import vals
from get_opt import *
from fit_opt_flares import calc_t90

# Retrieve the flare data
tab = get_full_opt()
isflare = tab['isflare'].values
nights_all = tab['mjdstart'].astype(int)
flare_nights_all = tab['mjdstart'][tab['isflare']==True].astype(int)
flare_nights = np.unique(tab['mjdstart'][tab['isflare']==True].astype(int)) 

# Iterate through all flares
for n,night in enumerate(flare_nights):
    dat = tab[nights_all==night]
    isflare_n = dat['isflare'].values

    # For each telescope that has a flare 
    flare_tels = np.unique(dat['#instrument'][isflare_n].values)
    for flare_tel in flare_tels:
        # Select all data from that telescope that night
        dat_tel = dat[dat['#instrument'].values==flare_tel]
        isflare_tel = dat_tel['isflare'].values

        # Get the time of the brightest detection
        flux = dat_tel['flux_extcorr'].values
        indpeak = np.argmax(flux[isflare_tel])
        tpeak = dat_tel['mjdstart'][isflare_tel].values[indpeak]

        # Already have the name of the telescope

        # Get the filters in which the flares were detected
        filt = np.unique(dat_tel[isflare_tel]['flt'].values)
        filtstr = ''.join(filt)

        # Get the T90 in the first of the listed filters
        use_filt = filt[0]
        wl = vals.leff[use_filt] # in AA
        freq = 3E18 / wl
        x = dat_tel['mjdstart'].values
        ey = dat_tel['unc_extcorr'].values
        t90_d,lpeak = calc_t90(freq, x, flux, ey)
        t90 = np.round(t90_d*24*60, 1)
        erad = np.trapz(flux, (x-x[0])*86400) * 1E-6 * 1E-23 * 4 * np.pi * vals.dL_cm**2 * freq
        erad_str = r"%s \times 10^{%s}"%(str(erad)[0:3], str(erad)[-2:])
        lpeak_str = r"%s \times 10^{%s}"%(str(lpeak)[0:3], str(lpeak)[-2:])
        printstr = "%s & %s & $%s$ & %s & $%s$ & $%s$" %(
                np.round(tpeak,4), tel, filtstr, t90, lpeak_str, erad_str)
        print(printstr)
