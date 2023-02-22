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
    filt_crit = lc['flt']==filt
    use_lc_filt = lc[filt_crit]

    # Either it's deep enough because there was a detection brighter
    # than the threshold 
    detected = np.logical_and(lc['mag']<threshold, lc['sig']>=3)

    # Or it's deep enough because there was no detection, and the limiting
    # magnitude is fainter than the threshold
    not_detected = np.logical_and.reduce((lc['maglim']>threshold, lc['sig']<3,
                                          np.abs(lc['maglim'])<99))
    thresh_crit = np.logical_or(detected, not_detected)
    choose = np.logical_and(thresh_crit, filt_crit)

    nrows = len(use_lc)
    print("%s rows" %nrows)
    if nrows>0:
        if nrows<10:
          print(use_lc)

        # Exposures where there was a flare this bright or brighter
        exp_on = use_lc[np.logical_and(
                        use_lc['mag']<=threshold, use_lc['isflare'])]
        time_on = sum(exp_on['exp'])
        n_exp_on = len(exp_on)

        # Calculate the duty cycle in several different ways.

        # The ratio of the time on to the time off (exposure times)
        print("calculate simple ratio of time on to total time")
        frac_time_on = time_on / sum(use_lc['exp'].values.astype(float))
        print(frac_time_on)

        # The binomial statistic---number of exposures
        print("binomial statistic with number of exposures")
        out = binomtest(n_exp_on, len(use_lc))
        val = out.statistic
        conf_low = out.proportion_ci().low
        conf_high = out.proportion_ci().high
        print(val, conf_low, conf_high)

        #print(
        #    filt, threshold, np.round(frac_time_on,2), np.round(val, 2), 
        #    np.round(conf_low, 3), np.round(conf_high, 2))

    else:
        print(filt, threshold, "no exposures")
