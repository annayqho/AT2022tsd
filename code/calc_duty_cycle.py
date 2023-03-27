"""
Calculate the duty cycle of the flares for different limiting magnitudes
"""

from numpy import random
import math
from scipy.stats import binomtest
from get_opt import *


def get_flaring_lc(threshold):
    """ Get the LC in the range where you think flaring is happening,
    for a certain threshold """
    # Get the optical photometry
    full_lc = get_full_opt()

    # Drop the columns you won't need
    full_lc = full_lc.drop(
            ['#instrument', 'flux', 'unc', 'sig', 'flare', 'emag', 
             'flux_extcorr', 'unc_extcorr', 'mag_extcorr', 'maglim_extcorr', 
             'nobs', 'istransient'], axis=1)

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

    # Convert all filters to g-band using the corrections I infer from 
    # f_nu ~ nu^{-1.6}
    # don't use w-band because it's a very wide filter...not clear how to do
    # the corrections. and anyway there are only a few points.
    lc = lc[lc.flt != 'w']

    lc.loc[lc['flt']=='r', 'mag'] = lc['mag'][lc['flt']=='r'] + 0.80 
    lc.loc[lc['flt']=='r', 'flt'] = 'g'

    lc.loc[lc['flt']=='i', 'mag'] = lc['mag'][lc['flt']=='i'] + 1.33
    lc.loc[lc['flt']=='i', 'flt'] = 'g'

    lc.loc[lc['flt']=='u', 'mag'] = lc['mag'][lc['flt']=='u'] - 0.81
    lc.loc[lc['flt']=='u', 'flt'] = 'g'

    # Drop the g-band filter column
    lc = lc.drop('flt', axis=1)

    # For any row where the mag is fainter than the maglim, set mag to 99
    lc.loc[lc['mag']>lc['maglim'], 'mag'] = 99.0

    # To be valid for this analysis, limiting magnitude of the image has to be
    # fainter than this threshold. But you have to convert it from 3-sigma
    # to 5-sigma, since you use a 5-sigma threshold to define flares.
    conversion_value = 2.5*np.log10(5)-2.5*np.log10(3) # positive number
    maglims = lc['maglim'] - conversion_value # make lims shallower
    deep_enough = maglims>threshold
    use_lc = lc[deep_enough]

    use_lc['maglim_5sig'] = maglims

    return use_lc


def get_favg(use_lc, threshold):
    """ For a given limiting magnitude, return 
    f_avg = T_on / T_tot 

    Parameters:
    -----------
    threshold: the limiting magnitude

    Returns:
    --------
    favg: the average duty cycle
    """
    nrows = len(use_lc)
    n_sensitive_exposures = nrows
    print("Nexp, ", n_sensitive_exposures)
    tot_time = sum(use_lc['exp'].values.astype(float))
    print("Texp, ", tot_time/60)
    if nrows>0:
        # Exposures where there was a flare this bright or brighter
        exp_on = use_lc[np.logical_and(
                        use_lc['mag']<=threshold, use_lc['isflare'])]
        time_on = sum(exp_on['exp'])
        n_exp_on = len(exp_on)

        # The ratio of the time on to the total time
        frac_time_on = time_on / tot_time
        print("Ton/Texp, ", frac_time_on)
        return frac_time_on
    else:
        print("not enough data for this threshold")
        return None


def get_duration(maglim):
    """ Get the minimum and maximum flare duration above 
    a certain limiting magnitude 

    I'm basically assuming half the observed minimum,
    to twice the observed maximum.

    returns
    -------
    T_min: minimum flare duration in days
    T_max = maximum flare duration in days
    """

    if maglim==22.5:
        T_min = 60
        T_max = 3*60*60
        return T_min/86400, T_max/86400
    elif maglim==24:
        T_min = 60
        T_max = 3*60*60
        return T_min/86400, T_max/86400
    elif maglim==21:
        T_min = 10*60
        T_max = 40*60
        return T_min/86400, T_max/86400
    else:
        print("invalid threshold")
        return None


def print_table_for_paper():
    """ Print a table with all the favg values """
    print("Threshold & $N_\mathrm{exp}$ & $T_\mathrm{exp}$ & $T_\mathrm{on}$/Texp")#& $N_\mathrm{exp,on}/N_\mathrm{exp}$")
    # Choose the filter and limiting magnitude
    for threshold in np.arange(19.5, 25, step=0.5):
        print("%s & %s & %s & %s \\\\" %(
            threshold,n_sensitive_exposures,int(tot_time/60),
            '{:.2f}'.format(frac_time_on)))#,'{:.2f}'.format(val),


if __name__=="__main__":
    ### Select your parameters

    ## Set the threshold
    thresh = 24

    # Get the the relevant exposures 
    lc = get_flaring_lc(thresh)

    # Get basic parameters
    favg = get_favg(lc, thresh) # duty cycle, avg
    T_min, T_max = get_duration(thresh)

    ## Set which duration to use
    #T = T_min
    T = T_max

    ## Set which range to use

    # for thresh = 21 mag, Tmin
    #avg_flare_rates = np.logspace(-1, 1) #result:[0.72, 7.54]->[0.005,0.05]

    # for thresh = 21 mag, Tmax
    #avg_flare_rates = np.logspace(-1.5,0.5)#result:[0.046,]->[0.001,]
    #avg_flare_rates = np.logspace(0.5, 2.5)#result:[,3.82]->[,0.1]

    # for thresh = 22.5 mag, Tmin
    #avg_flare_rates = np.logspace(1,3)  # duty cycle: [0.038,]
    #avg_flare_rates = np.logspace(2,4)  # duty cycle: [,0.28]
    # for thresh = 22.5 mag, Tmax
    # avg_flare_rates = np.logspace(-2,0) # duty cycle: [0.012,]
    # avg_flare_rates = np.logspace(0,2)  # duty cycle: [,0.56]

    # for thresh = 24 mag, Tmin ... didn't run
    # because it seems to always depend on the long ones anyway
    # avg_flare_rates = np.logspace(-2,0)  # 
    # avg_flare_rates = np.logspace(0,2)  # 

    # for thresh = 24 mag, Tmax 
    #avg_flare_rates = np.logspace(-1,1)  # 0.029
    avg_flare_rates = np.logspace(0,2)  #  max is > 1...

    # To convert from a flare rate to a duty cycle, multiply by the 
    # fraction of a day taken up by one flare. 

    ### Don't change anything below

    # Construct a set of burst times that obey a Poisson distribution
    # For Poisson, the time between events is exponentially distributed
    tstart = min(lc['mjdstart'])-1 # one day before the start of the window
    tend = max(lc['mjdstart'])+1  # one day after the end of the window

    for j,avg_flare_rate in enumerate(avg_flare_rates):
        print("running for %s mag and duration %s min" %(thresh, T*24*60))
        print("rate %s" %avg_flare_rate)
        print("duty cycle %s" %str(avg_flare_rate*T))
        duty_cycles = []
        for i in np.arange(1000):
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
            else:
                duty_cycles.append(0)
        duty_cycles = np.array(duty_cycles)
        frac_high = sum(duty_cycles>favg)/len(duty_cycles) 
        frac_low = sum(duty_cycles<favg)/len(duty_cycles) 
        if np.logical_or(frac_high>0.975, frac_low>0.975):
            print("not consistent")
        else:
            print("consistent")
