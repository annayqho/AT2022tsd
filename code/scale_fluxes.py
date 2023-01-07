""" Scale fluxes using variations in 1635 at that frequency

Notes that I deleted Jul 23 from the file and June 28
and deleted comments from July 2
"""

import numpy as np
from astropy.time import Time

def sma_lc():
    data_dir = "/Users/annaho/Dropbox/astro/papers/papers_complete/AT2018cow/data"
    lines = open("%s/SMA_AT2018cow_quasar.txt" %data_dir, "r").readlines()
    t0 = Time(58285, format='mjd') # our zero-point
    dt = []
    nu = []
    f = []
    ef = []
    ref = []

    for line in lines:
        has_date = 'Jun' in line or 'Jul' in line or 'Aug' in line
        items = line.split()
        if has_date and len(items) > 10:
            date = items[1]
            if date == 'Jun':
                mm = '06'
            elif date == 'Jul':
                mm = '07'
            elif date == 'Aug':
                mm = '08'
            else:
                print("something went wrong")

            t_raw = Time('2018-%s-%s' %(mm,items[2])).mjd + float(items[3])/24
            dt.append(t_raw - t0.mjd)
            nu.append(line.split()[6])
            f.append(line.split()[7])
            ef.append(line.split()[9])
            ref.append(line.split()[10])
    dt = np.array(dt, dtype=float)
    nu = np.array(nu, dtype=float)
    f = np.array(f, dtype=float)
    ef = np.array(ef, dtype=float)
    ref = np.array(ref, dtype=float)

    # Step 1: Measure the average quasar flux at each frequency
    u_nu = np.unique(nu)
    u_ref = np.array([np.mean(ref[nu==val]) for val in u_nu])

    # Step 2: Scale each flux measurement by the 
    # ratio of this avg/1635 flux at that epoch
    npts = len(dt)
    f_scaled = np.zeros(npts)
    ef_scaled = np.zeros(npts)

    for ii,fmeas in enumerate(f):
        mean_ref = u_ref[u_nu==nu[ii]][0]
        ratio = mean_ref / ref[ii]
        f_scaled[ii] = f[ii] * ratio
        ef_scaled[ii] = ef[ii] * ratio

    # Generate the light curve in the range 230.6-234.6
    dt_f_230 = []
    flux_f_230 = []
    eflux_f_230 = []

    # Count the number of points on this day measured at 231.5 GHz
    for ii,dt_val in enumerate(np.unique(dt)):
        choose = np.logical_and(
                nu[dt==dt_val] <= 234.6, nu[dt==dt_val] >= 230.6)
        if sum(choose) > 0:
            # then add this point
            dt_f_230.append(dt_val)
        # If the answer is 1, then store that point
        if sum(choose) == 1:
            flux_f_230.append(f_scaled[dt==dt_val][choose][0])
            eflux_f_230.append(ef_scaled[dt==dt_val][choose][0])
        # If the answer is 2, then take the average of those points
        if sum(choose) == 2:
            w = 1/(ef_scaled[dt==dt_val][choose])**2
            fmean, wsum = np.average(
                    f_scaled[dt==dt_val][choose], 
                    weights=w,
                    returned=True)
            efmean = np.sqrt(1/np.sum(w))
            flux_f_230.append(fmean)
            eflux_f_230.append(efmean)
        # If the answer is > 2, something is wrong...
        if sum(choose) > 2:
            print("something is wrong!")
        if sum(choose) == 0:
            # the only days that aren't within this tiny range
            # have nu = 231.5 and nu = 218.0,
            # so these are the ones I will bother to scale
            # choose the points with 243.3 and scale them to 231.5
            toscale = nu[dt==dt_val] == 243.3
            if sum(toscale) == 1:
                dt_f_230.append(dt_val)
                flux_f_230.append(f_scaled[dt==dt_val][toscale][0] * (243.3/231.5))
                eflux_f_230.append(
                        ef_scaled[dt==dt_val][toscale][0] * (243.3/231.5))
            else:
                toscale = nu[dt==dt_val] == 218
                if sum(toscale) == 1:
                    dt_f_230.append(dt_val)
                    flux_f_230.append(f_scaled[dt==dt_val][toscale][0] * (218.0/231.5))
                    eflux_f_230.append(
                            ef_scaled[dt==dt_val][toscale][0] * (218.0/231.5))
    dt_f_230 = np.array(dt_f_230)
    flux_f_230 = np.array(flux_f_230)
    eflux_f_230 = np.array(eflux_f_230)

    
    # Generate the light curve in the range 341.5-349 GHz
    dt_f_345 = []
    flux_f_345 = []
    eflux_f_345 = []
    for ii,dt_val in enumerate(np.unique(dt)):
        # Count the number of points on this day measured near 345 GHz
        choose = np.logical_and(
                nu[dt==dt_val] >= 341.5, nu[dt==dt_val] <= 349)
        if sum(choose) > 0:
            # then add this point
            dt_f_345.append(dt_val)
        # If the answer is 1, then store that point
        if sum(choose) == 1:
            flux_f_345.append(f_scaled[dt==dt_val][choose][0])
            eflux_f_345.append(ef_scaled[dt==dt_val][choose][0])
        # If the answer is 2, then take the average of those points
        if sum(choose) == 2:
            w = 1/(ef_scaled[dt==dt_val][choose])**2
            fmean, wsum = np.average(
                    f_scaled[dt==dt_val][choose], 
                    weights=w,
                    returned=True)
            efmean = np.sqrt(1/np.sum(w))
            flux_f_345.append(fmean)
            eflux_f_345.append(efmean)
        # If the answer is > 2, something is wrong...
        if sum(choose) > 2:
            print("something is wrong!")
    dt_f_345 = np.array(dt_f_345)
    flux_f_345 = np.array(flux_f_345)
    eflux_f_345 = np.array(eflux_f_345)

    return (dt,nu,f_scaled,ef_scaled), (dt_f_230, flux_f_230, eflux_f_230), (dt_f_345, flux_f_345, eflux_f_345)
