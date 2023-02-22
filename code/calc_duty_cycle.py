"""
Calculate the duty cycle of the flares for different limiting magnitudes
"""

from scipy.stats import binomtest
from get_opt import *

# Get the optical photometry
full_lc = get_full_opt()

# I think it's convincing that the flaring starts with the first
# ZTF detection (aka the first flare detection).
min_mjd = np.min(full_lc['mjdstart'][full_lc['isflare']])

# Let's assume the flaring ends with the last flare detection.
max_mjd = np.max(full_lc['mjdstart'][full_lc['isflare']])

# Use the LC in the relevant range
lc = full_lc[np.logical_and(
    full_lc['mjdstart']>=min_mjd, full_lc['mjdstart']<=max_mjd)]
lc = lc.sort_values(by=['mjdstart'], ignore_index=True, axis=0)

# Let R = r
lc['flt'][lc['flt']=='R'] = ['r']*len(lc[lc['flt']=='R'])

for filt in np.unique(lc['flt']):
    print("Running for %s" %filt)
    # Choose the filter and limiting magnitude
    threshold = 21

    # For a given filter and limiting magnitude, 
    # identify all sufficiently deep exposures
    thresh_crit = np.logical_or(
            np.logical_and(lc['mag']>threshold, lc['mag']<99), 
            np.logical_and(lc['maglim']>threshold, lc['maglim']<99))
    filt_crit = lc['flt']==filt
    choose = np.logical_and(thresh_crit, filt_crit)
    lc = lc[choose]

    print("%s rows" %len(lc))
    if len(lc)>0:
        # Calculate the duty cycle in several different ways.

        # The ratio of the time on to the time off (exposure times)
        frac_time_on = sum(lc['exp'][lc['isflare']==True])/sum(lc['exp'])

        # The binomial statistic---number of exposures
        out = binomtest(len(lc['exp'][lc['isflare']==True]),len(lc['exp']))
        val = out.statistic
        conf_low = binomtest(
            len(lc['exp'][lc['isflare']==True]),len(lc['exp'])).proportion_ci().low
        conf_high= binomtest(
            len(lc['exp'][lc['isflare']==True]),len(lc['exp'])).proportion_ci().high

        print(
            filt, threshold, np.round(frac_time_on,2), np.round(val, 2), 
            np.round(conf_low, 3), np.round(conf_high, 2))

    else:
        print(filt, threshold, "no exposures")
