#from astropy.timeseries import LombScargle
import matplotlib.pyplot as plt
import cmasher as cmr
from get_opt import *

dat = get_full_opt()
#choose = np.logical_and.reduce((np.logical_or.reduce((dat['flt']=='r', dat['flt']=='g', dat['flt']=='i', dat['flt']=='w')), dat['sig']>-99, dat['mjdstart']>59856.4))

choose = np.logical_and.reduce((np.logical_or.reduce((dat['flt']=='r', dat['flt']=='g', dat['flt']=='i', dat['flt']=='w')), dat['sig']>-99, dat['mjdstart']>59856.4, dat['mjdstart']<59943.0))

t = dat['mjdstart'][choose].values.astype(float)
f = dat['flux'][choose].values.astype(float)
ef = dat['unc'][choose].values.astype(float)
isf = dat['isflare'][choose].values
filt = dat['flt'][choose].values
tel = dat['#instrument'].values[choose]

#ls = LombScargle(t, f, ef)
#frequency, power = ls.autopower(method='slow', samples_per_peak=5)
#period = 1/frequency
#plot(period, power, 'k')

#period_sorted = period[np.argsort(power)][::-1]
period_sorted = np.arange(3/24, 1, 5/60/24)

cols = cmr.take_cmap_colors(
    'cmr.rainforest', 11, cmap_range=(0.1, 0.9), return_fmt='hex')[::-1]
m = ['o', 'D', 's', '*', 'v', '>', '<', '^', 'H', 'o', 'D']

for i in np.arange(len(period_sorted)):
    fig,ax = plt.subplots(1,1,figsize=(8,4))
    P = period_sorted[i]
    print(P)

    for j,tel_val in enumerate(np.unique(tel)):
        choose = np.logical_and(tel==tel_val, isf)
        plt.errorbar(t[choose]/P % 1, f[choose], ef[choose], fmt=m[j], 
                     ms=4,label=tel_val,
                     c=cols[j], zorder=2)
        #choose = np.logical_and(tel==tel_val, ~isf)
        #plt.errorbar(t[choose]/P % 1, f[choose], ef[choose], fmt=m[j], ms=4,alpha=0.1,
        #             c=cols[j], zorder=0)
    plt.title(P*24)
    plt.legend()
    plt.xlim(-0.1,1.1)
    plt.savefig("period_search_3hr_1d_5minstep/P_only_flare_period_%s.png" %(i))
    plt.close()

    # Tentative period at 5.88 hours?
