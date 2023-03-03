"""
Calculate the duty cycle of the flares for different limiting magnitudes
"""

from numpy import random
import math
from scipy.stats import binomtest
from get_opt import *


def get_flaring_lc(filt, threshold):
    """ Get the LC in the range where you think flaring is happening,
    for a certain filter and threshold """
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

    # Update name of filter
    lc['flt'][lc['flt']=='R'] = ['r']*len(lc[lc['flt']=='R'])

    # For the filter and limiting magnitude, identify all relevant exposures
    filt_crit = lc['flt']==filt
    use_lc_filt = lc[filt_crit]

    # To be valid for this analysis, limiting magnitude of the image has to be
    # fainter than this threshold. But you have to convert it from 3-sigma
    # to 5-sigma, since you use a 5-sigma threshold to define flares.
    conversion_value = 2.5*np.log10(5)-2.5*np.log10(3) # positive number
    maglims = use_lc_filt['maglim'] - conversion_value # make lims shallower
    deep_enough = maglims>threshold
    use_lc = use_lc_filt[deep_enough]

    return use_lc


def get_favg(use_lc, filt, threshold):
    """ For a given filter and limiting magnitude, return 
    f_avg = T_on / T_tot 

    Parameters:
    -----------
    filt: filter ('r', 'g', etc)
    threshold: the limiting magnitude

    Returns:
    --------
    favg: the average duty cycle
    """
    nrows = len(use_lc)
    n_sensitive_exposures = nrows
    tot_time = sum(use_lc['exp'].values.astype(float))
    if nrows>0:
        # Exposures where there was a flare this bright or brighter
        exp_on = use_lc[np.logical_and(
                        use_lc['mag']<=threshold, use_lc['isflare'])]
        time_on = sum(exp_on['exp'])
        n_exp_on = len(exp_on)

        # The ratio of the time on to the total time
        frac_time_on = time_on / tot_time
        return frac_time_on
    else:
        print("not enough data for this filter/threshold")
        return None


def get_duration():
    """ Get the minimum and maximum flare duration above 
    a certain limiting magnitude 

    returns
    -------
    T_min: minimum flare duration in days
    T_max = maximum flare duration in days
    """

    # For now, let's just adopt some typical minimum and maximum
    # value, since the diversity of morphologies is highly uncertain.
    # The OBSERVED minimum is a couple of minutes, so let's say 1 min.
    # The OBSERVED maximum is a copule of hours, so let's say 3 hours.

    T_min = 60
    T_max = 3*60*60
    return T_min/86400, T_max/86400


def print_table_for_paper():
    """ Print a table with all the favg values """
    print("Band & Threshold & $N_\mathrm{exp}$ & $T_\mathrm{exp}$ & $T_\mathrm{on}$/Texp")#& $N_\mathrm{exp,on}/N_\mathrm{exp}$")
    for filt in np.unique(lc['flt']):
        #print("Running for %s" %filt)
        # Choose the filter and limiting magnitude
        for threshold in np.arange(19.5, 25, step=0.5):
            print("$%s$ & %s & %s & %s & %s \\\\" %(
                filt,threshold,n_sensitive_exposures,int(tot_time/60),
                '{:.2f}'.format(frac_time_on)))#,'{:.2f}'.format(val),


if __name__=="__main__":
    ### Select your parameters
    filt = 'g'
    thresh = 22

    # Get the the relevant exposures 
    lc = get_flaring_lc(filt,thresh)

    # Get basic parameters
    favg = get_favg(lc, filt, thresh) # duty cycle, avg
    T_min, T_max = get_duration()

    ### Which duration to use
    T = T_max # bounds: 0.2/day-3.6/day, duty cycle 0.025-0.45
    avg_flare_rates = np.logspace(-1,1)

    T = T_min # bounds: 87-290/day, duty cycle 0.06-0.20
    #avg_flare_rates = np.logspace(1,2)

    # Construct a set of burst times that obey a Poisson distribution
    # For Poisson, the time between events is exponentially distributed
    tstart = min(lc['mjdstart'])-1 # one day before the start of the window
    tend = max(lc['mjdstart'])+1  # one day after the end of the window

    for j,avg_flare_rate in enumerate(avg_flare_rates):
        print("rate %s" %avg_flare_rate)
        duty_cycles = []
        for i in np.arange(600):
            print(i)
            # Arrival times
            r = random.rand(100000)
            inter_flare_times = -np.log(r)/avg_flare_rate

            # Flare start times
            flare_start_times = tstart+np.cumsum(inter_flare_times)
            flare_start_times = flare_start_times[flare_start_times<tend]

            if len(flare_start_times)>0:
                # Flare end times, given the duration
                flare_end_times = flare_start_times + T

                # Calculate the duty cycle given our observing times
                obs_times = lc['mjdstart'].values
                obs_exp = lc['exp'].values
                ton = 0
                for i,t in enumerate(obs_times):
                    det = np.logical_and(t>=flare_start_times, t<=flare_end_times)
                    if sum(det)>0:
                        ton += obs_exp[i]
                duty_cycle = ton / sum(obs_exp)
                duty_cycles.append(duty_cycle)
        duty_cycles = np.array(duty_cycles)
        frac_high = sum(duty_cycles>favg)/len(duty_cycles) 
        frac_low = sum(duty_cycles<favg)/len(duty_cycles) 
        if np.logical_or(frac_high>0.975, frac_low>0.975):
            print("not consistent")
        else:
            print("consistent")
