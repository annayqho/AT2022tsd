""" Calculate the properties of each flare: duration, luminosity, etc """

import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
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

        # Sometimes, there's just one flare detection.
        if len(dat_tel[isflare_tel])==1:
            tpeak = dat_tel[isflare_tel]['mjdstart'].values[0]
            tel = dat_tel[isflare_tel]['#instrument'].values[0]
            fpeak = dat_tel[isflare_tel]['flux_extcorr'].values[0]
            filtstr = dat_tel[isflare_tel]['flt'].values[0]
            wl = vals.leff[filtstr] # in AA
            freq = 3E18 / wl
            lpeak = fpeak * 1E-6 * 1E-23 * 4 * np.pi * vals.dL_cm**2 * freq
            t90 = '--' # can't calculate
            erad_str = '--'
        else:
            flux = dat['flux_extcorr'].values
            indpeak = np.argmax(flux)
            tpeak = dat['mjdstart'].values[indpeak]
            tel = np.unique(dat['#instrument'].values)[0]
            filt = np.unique(dat['flt'].values)
            filtstr = ''.join(filt)
            use_filt = filt[0]
            wl = vals.leff[use_filt] # in AA
            freq = 3E18 / wl
            t90_d,lpeak = calc_t90(freq, dat['mjdstart'].values,
                                 flux, dat['unc_extcorr'].values)
            t90 = np.round(t90_s*24*60, 1)
            erad = t90*60*lpeak
            erad_str = r"%s \times 10^{%s}"%(str(erad)[0:3], str(erad)[-2:])
        lpeak_str = r"%s \times 10^{%s}"%(str(lpeak)[0:3], str(lpeak)[-2:])
        printstr = "%s & %s & $%s$ & %s & $%s$ & $%s$" %(
                np.round(tpeak,4), tel, filtstr, t90, lpeak_str, erad_str)
        print(printstr)
