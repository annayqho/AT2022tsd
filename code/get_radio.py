""" Get the radio data from SMA, ATCA, ALMA """

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from synchrotron_fit import fit_self_abs, self_abs


def get_data_all():
    """ 
    Get the full set of data: 

    Returns
    -------
    tel: array telling you which telescope the data point is from
    freq: freq of the data point
    days: day of the data point, from time of explosion
    flux: flux of the data point
    eflux_form: formal uncertainty in flux of the data point
    eflux_sys: estimate of systematic uncertainty in flux of the data point
    """
    data_dir = "../data"
    dat = Table.read(
        "%s/radio_lc.dat" %data_dir, delimiter="&", format='ascii.no_header')
    days = np.array(dat['col1'])
    tel = np.array(dat['col2'])
    freq = np.array(dat['col3']).astype(float)
    flux_raw = np.array(dat['col4'])
    flux = np.zeros(len(flux_raw))
    eflux_form = np.zeros(len(flux_raw))
    eflux_sys = np.zeros(len(flux_raw))

    # Now, go through and get the appropriate uncertainties
    for ii,val in enumerate(flux_raw):
        if tel[ii] == 'ALMA':
            # ALMA points only have systematic uncertainties
            add = val.split("$")[1].split('pm')[0]
            flux[ii] = add
            eflux_form[ii] = -99
            if freq[ii] < 200: # bands 3 and 4
                eflux_sys[ii] = 0.05*flux[ii]
            elif freq[ii] > 600: # band 9
                eflux_sys[ii] = 0.20*flux[ii]
            else: # bands 7 and 8
                eflux_sys[ii] = 0.10*flux[ii]
        elif tel[ii] == 'SMA':
            # SMA points have formal uncertainties
            flux[ii] = float(val.split("pm")[0][1:])
            eflux_form[ii] = float(val.split("pm")[1][0:-1])
            # 10% sysetmatic uncertainty
            eflux_sys[ii] = 0.10*flux[ii]
        elif tel[ii] == 'ATCA':
            # ATCA points have both formal and systematic uncertainties
            flux[ii] = float(val.split("pm")[0][1:])
            eflux_form[ii] = float(val.split("pm")[1][0:-1])
            eflux_sys[ii] = 0.10*flux[ii]
    return tel, freq, days, flux, eflux_form, eflux_sys 


def get_spectrum(day):
    """ 
    Get full for a specific day
    For each frequency band, you will have to interpolate.
    You're also allowed to extrapolate if you're at lower frequencies,
    since the self-absorbed parts are well-behaved.

    For the SMA data, scale the data according to variation
    in the reference quasar.
    """
    # First, get *all* of the data
    tel, freq, days, flux, eflux_form, eflux_sys = get_data_all()
    print(day)
    # Don't generate eflux_tot here,
    # because the right thing to do depends on the telescope

    # The way we interpolate depends on the telescope.
    # For ATCA data, you can assume that the data is rising as a power law
    # in log-log space.
    # The frequencies are 5.5, 9, 16.7, 21.2, and 34 GHz.
    # ignore upper limits
    choose_tel = np.logical_and(tel == 'ATCA', eflux_form > 0)
    ufreq = np.unique(freq[choose_tel])
    nus = [] # which frequencies to keep
    ms = [] # power-law fit
    bs = [] # power-law fit
    
    for val in ufreq:
        # fit a power law in log-log space
        choose = np.logical_and(choose_tel, freq==val)
        # no upper limits

        # only bother if there's >1 point at this frequency
        if sum(choose) > 1:
            nus.append(val)
            # if you allow weights, then the 34 GHz points look really weird
            m,b = np.polyfit(
                    np.log10(days[choose]), np.log10(flux[choose]), deg=1)#,
                    #w=1/eflux_tot[choose]**2)
            ms.append(m)
            bs.append(b)
    
    # Now, interpolate for that day to get the spectrum
    spec = []
    for ii,val in enumerate(nus):
        # Get light curve for that day
        spec.append(10**(ms[ii] * np.log10(day) + bs[ii]))

    # OK, "nus" and "spec" have been updated for the ATCA band!

    # Next: SMA.
    # This is much more complicated, since the emission is not self-absorbed
    # and the frequencies are not consistent.
    
    # I think the right thing to do is to interpolate the data
    # for each frequency...

    # The frequencies with by far the most data are 215.5 GHz and 231.5 GHz
    # And the one with definitively the most data is 231.5

    # So I think the right thing to do is show the 231.5 GHz light curve.
    # For every day that wasn't observed with the 231.5 GHz receiver,
    # interpolate the spectrum if you can

    # the higher frequencies are more difficult.
    # on each day, there are at least two measurements
    # that span 345 GHz.
    # so, on each day interpolate to estimate the flux at 345 GHz.
    tel, freq, days, flux, eflux_form, eflux_sys = get_data_all()
    eflux_tot = np.sqrt(eflux_form**2 + eflux_sys**2)

    choose_tel = tel == 'SMA' 
    d = [] # day
    f_lo = [] # interpolated flux at 231.5 GHz
    f_hi = [] # interpolated flux at 345 GHz
    uday = np.unique(days)

    for ii,val in enumerate(uday):
        choose = np.logical_and(choose_tel, days==val)
        success = False # at least one point added

        # lower bands
        keep_flux = flux[choose][freq[choose] == 231.5]
        if len(keep_flux) == 1:
            # Then you can save the value directly
            success = True
            f_lo.append(keep_flux[0])
        elif len(keep_flux) > 1:
            print("more than one point on a given day?")

        else:
            # Then you can try to interpolate
            try:
                xfit = freq[choose]
                yfit = flux[choose]
                order = np.argsort(xfit)
                out = np.interp(231.5, xfit[order], yfit[order])
                f_lo.append(out)
                success = True
            except:
                # don't extrapolate
                pass

        # upper bands
        keep_flux = flux[choose][freq[choose] == 345]
        if len(keep_flux) == 1:
            # Then you can save the value directly
            f_hi.append(keep_flux[0])
            success = True
        else:
            # Then you can try to interpolate
            try:
                xfit = freq[choose]
                yfit = flux[choose]
                order = np.argsort(xfit)
                out = np.interp(345, xfit[order], yfit[order])
                f_hi.append(out)
                success = True
            except:
                # don't extrapolate
                pass
        if success:
            d.append(val)

    # Now, interpolate for that day to get the two-point spectrum
    nus.append(231.5)
    d = np.array(d)
    f_lo = np.array(f_lo)
    f_hi = np.array(f_hi)
    order = np.argsort(d)
    spec.append(np.interp(day, d[order], f_lo[order]))
    nus.append(345)
    spec.append(np.interp(day, d[order], f_hi[order]))

    # Finally: ALMA
    # For any ALMA points, you should just return the value on the day
    # if it exists.
    # Nothing to interpolate.
    choose = np.logical_and(days == day, tel == 'ALMA')
    for ii,nuval in enumerate(freq[choose]):
        nus.append(nuval)
        spec.append(flux[choose][ii])

    return np.array(nus), np.array(spec)
