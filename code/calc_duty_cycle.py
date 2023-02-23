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

print("Band & Threshold & $N_\mathrm{exp}$ & $T_\mathrm{exp}$ & $T_\mathrm{on}$/Texp")#& $N_\mathrm{exp,on}/N_\mathrm{exp}$")
for filt in np.unique(lc['flt']):
    #print("Running for %s" %filt)
    # Choose the filter and limiting magnitude
    for threshold in np.arange(19.5, 25, step=0.5):
        # For a given filter and limiting magnitude, 
        # identify all sufficiently deep exposures
        filt_crit = lc['flt']==filt
        use_lc_filt = lc[filt_crit]

        # To be valid for this analysis, EITHER...

        # (a) the limiting magnitude of the image has to be 
        # fainter than this threshold. But you have to convert it from 3-sigma
        # to 5-sigma, since you use a 5-sigma threshold for your flares.
        conversion_value = 2.5*np.log10(5)-2.5*np.log10(3)
        maglims = use_lc_filt['maglim'] - conversion_value
        deep_enough = maglims>threshold

        # OR (b) there's a flare detection that's bright enough (or both)
        # no, I don't think that's correct. you should only use the images
        # that are sensitive enough, because otherwise for very faint magnitudes
        # you end up way overestimating the flare duty cycle.
        #bright_enough = np.logical_and(
        #        use_lc_filt['mag']<threshold, use_lc_filt['isflare'])

        #use_lc = use_lc_filt[np.logical_or(deep_enough, bright_enough)]
        use_lc = use_lc_filt[deep_enough]

        nrows = len(use_lc)
        n_sensitive_exposures = nrows
        tot_time = sum(use_lc['exp'].values.astype(float))
        #print("%s rows" %nrows)
        if nrows>0:
            #if nrows<10:
              #print(use_lc)

            # Exposures where there was a flare this bright or brighter
            exp_on = use_lc[np.logical_and(
                            use_lc['mag']<=threshold, use_lc['isflare'])]
            time_on = sum(exp_on['exp'])
            n_exp_on = len(exp_on)

            # Calculate the duty cycle in several different ways.

            # The ratio of the time on to the total time
            #print("calculate simple ratio of time on to total time")
            frac_time_on = time_on / tot_time
            #print(frac_time_on)

            # The binomial statistic---number of exposures
            #print("binomial statistic with number of exposures")
            out = binomtest(n_exp_on, nrows)
            val = out.statistic
            conf_low = out.proportion_ci().low
            conf_high = out.proportion_ci().high
            #print(val, conf_low, conf_high)

            #print("%s & %s & %s & %s & %s & %s [%s, %s]" %(
            #    filt,threshold,n_sensitive_exposures,int(tot_time/60),
            #    '{:.2f}'.format(frac_time_on),'{:.2f}'.format(val),
            #    '{:.2f}'.format(conf_low),'{:.2f}'.format(conf_high)))
            print("$%s$ & %s & %s & %s & %s \\\\" %(
                filt,threshold,n_sensitive_exposures,int(tot_time/60),
                '{:.2f}'.format(frac_time_on)))#,'{:.2f}'.format(val),
                #'{:.2f}'.format(conf_low),'{:.2f}'.format(conf_high)))

        #else:
            #print(filt, threshold, "no exposures")
