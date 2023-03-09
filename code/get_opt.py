""" Code to get the optical photometry """

import pandas as pd
import numpy as np
import vals
from astropy.time import Time


ddir = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/opt"


def measure_baseline(fcqf_use, ref_start, ref_end):
    """ Measure the baseline flux level for a given field-ccd-quadrant-filter
    combination that uses the same reference image """
    dat = pd.read_table(
            "%s/ipac_forced_phot_lc_full_survey.txt" %ddir, comment='#', 
            delimiter=' ')
    fid = np.zeros(len(dat['filter,']))
    fid[dat['filter,']=='ZTF_g'] = 1
    fid[dat['filter,']=='ZTF_r'] = 2
    fid[dat['filter,']=='ZTF_i'] = 3
    fcqf = (dat['field,']*10000+dat['ccdid,']*100+dat['qid,']*10+fid).astype(int).values
    choose = np.logical_and.reduce(
            (fcqf==fcqf_use, dat['refjdend,']==ref_end, 
             dat['refjdstart,']==ref_start))
    #x = dat['jd,'].values[choose]
    y = dat['forcediffimflux,'].values[choose]
    #ey = dat['forcediffimfluxunc,'].values[choose]
    #plt.errorbar(x[x<t0-10],y[x<t0-10],ey[x<t0-10],fmt='o')
    return np.median(y[~np.isnan(y)])


def get_ipac(inputf="%s/ipac_forced_phot.txt" %ddir):
    """ Get IPAC forced photometry """
    # Extract table
    a = pd.read_table(inputf, comment='#', delimiter=' ')

    # Extract all field/ccd/filter combinations
    fid = np.zeros(len(a['filter,']))
    fid[a['filter,']=='ZTF_g'] = 1
    fid[a['filter,']=='ZTF_r'] = 2
    fid[a['filter,']=='ZTF_i'] = 3
    fcqf = (a['field,']*10000+a['ccdid,']*100+a['qid,']*10+fid).astype(int).values
    ufcqf = np.unique(fcqf) # 5060333, 5060332, 5060331

    # Reference for i-band ends 30 Sept 2022...but I think it's OK,
    # since the references are median stacks.
    #refends = np.unique(a['refjdend,'][fcqf==5060333].values)
    #iband = fid==3
    #a = a[~iband]
    #fcqf = fcqf[~iband]

    # Bad pixel contamination: 4 images
    if a['procstatus'].dtype=='int64':
        bad_pix = a['procstatus'].values == 56
    else:
        bad_pix = np.logical_or(a['procstatus'].values == '56',
            a['procstatus'].values == '56,57')

    # Remove two of the images since you can't inspect them manually on IPAC
    # (maybe not publicly available yet)
    bad_jd = np.logical_and(bad_pix, a['jd,']>2459840)
    a = a[~bad_jd]
    fcqf = fcqf[~bad_jd]

    # Bad observing conditions
    # Criteria identified by plotting zpmaginpsci, zpmaginpscirms, and scisigpix
    # against each other
    bad_conditions = a['zpmaginpsci,'].values < 25.5
    a = a[~bad_conditions]
    fcqf = fcqf[~bad_conditions]

    # Measure the baseline flux level for each combo
    baseline_g = measure_baseline(5060331, 2458323.943912, 2458431.907211)
    baseline_r = measure_baseline(5060332, 2458340.956678, 2458437.829097)

    # Get the full light curve
    jd = a['jd,'].values
    flux = a['forcediffimflux,'].values
    eflux = a['forcediffimfluxunc,'].values
    filt = a['filter,'].values
    maglim = a['diffmaglim,'].values
    zp = a['zpdiff,'].values
    exp = a['exptime,'].values

    # Update the name of the filters
    filt[filt=='ZTF_g'] = 'g'
    filt[filt=='ZTF_r'] = 'r'
    filt[filt=='ZTF_i'] = 'i'

    # Subtract the baseline values
    flux[fcqf==5060331] = flux[fcqf==5060331]-baseline_g
    flux[fcqf==5060332] = flux[fcqf==5060332]-baseline_r

    # flux values are NaN? NOPE
    # nan = np.isnan(flux)

    # some chisq values are bad
    # NONE
    # bad_chisq = np.isnan(a['forcediffimchisq,'])

    # Look at the forcediffimchisq and make sure the average is 1
    # and that it's not dependent on flux
    forcediffimchisq = a['forcediffimchisq,'].values

    # Cobble together the LC
    SNT = 3 # threshold for declaring a measurement a non-detection
    SNU = 3 # S/N to use when computing an upper limit

    mag = np.array([99]*len(flux)).astype(float)
    emag = np.array([99]*len(eflux)).astype(float)
    is_det = flux/eflux >= SNT
    mag[is_det] = zp[is_det]-2.5*np.log10(flux[is_det])
    emag[is_det] = 1.0857*eflux[is_det]/flux[is_det]

    # Correct for MW extinction
    mag_extcorr = np.array([99]*len(mag)).astype(float)
    for f in np.unique(filt):
        choose = np.logical_and(is_det, filt==f)
        mag_extcorr[choose] = mag[choose]-vals.ext[f]

    # Convert the flux and eflux into physical units, uJy
    f0 = 10**(0.4*zp)
    fratio = flux/f0
    fujy = fratio * 3631 * 1E6
    efujy = eflux * (fujy/flux)

    # Correct for MW extinction
    fujy_extcorr = np.ones(len(fujy)).astype(float)
    efujy_extcorr = np.ones(len(fujy)).astype(float)
    for f in np.unique(filt):
        choose = filt==f
        fac = 10**(vals.ext[f]/2.5)
        fujy_extcorr[choose] = fujy[choose] * fac
        efujy_extcorr[choose] = efujy[choose] * fac

    # Final arrays for the light curve
    order = np.argsort(jd)
    jd = jd[order]
    fujy = fujy[order]
    efujy = efujy[order]
    fujy_extcorr = fujy_extcorr[order]
    efujy_extcorr = efujy_extcorr[order]
    mag = mag[order]
    mag_extcorr = mag_extcorr[order]
    emag = emag[order]
    exp = exp[order]
    filt = filt[order]

    return jd,exp,filt,mag,mag_extcorr,emag,\
           fujy,efujy,fujy_extcorr,efujy_extcorr


def get_last():
    """ Get observations from the LAST telescope. Data table provided
    by Eran """
    dat = pd.read_fwf(ddir+"/Table_LAST_binned.txt")
    mjd = Time(dat['% JD-2459000']+2459000, format='jd').mjd
    flux = dat['Flux,'].values
    eflux = dat['FluxErr,'].values
    nobs = dat['Nobs'].values
    return mjd,flux,eflux,nobs


def get_lulin():
    """ Get Lulin observations. Limit is 3-sigma. 
    """
    dat = pd.read_csv(ddir + "/AT2022tsd_LOT_SLT.txt", delimiter=' ',
                      names=['MJD','maglim','filt'])
    return dat


def get_panstarrs():
    """ Get the Pan-STARRS data """
    inputf = ddir + "/AT2022tsd-PS1-forced.csv"
    dat  = pd.read_csv(inputf)

    # Add a magnitudes column
    dat['mag'] = [99]*len(dat)
    dat['emag'] = [99]*len(dat)
    dat['maglim'] = [99]*len(dat)

    # Detections
    SNU = 3 # provide 3-sigma U.L.
    SNT = 3
    dat['sig'] = dat['ujy'].values/dat['dujy'].values
    isdet = dat['sig']>=SNT # confident detection
    fdet = dat['ujy'][isdet]
    efdet = dat['dujy'][isdet]

    # Correct for MW extinction
    f_extcorr = np.copy(dat['ujy'].values)
    ef_extcorr = np.copy(dat['dujy'].values)
    for f in np.unique(dat['filter']):
        choose = dat['filter'].values==f
        fac = 10**(vals.get_extinction(vals.ps1_leff[f])/2.5)
        f_extcorr[choose] = dat['ujy'].values[choose]*fac
        ef_extcorr[choose] = dat['dujy'].values[choose]*fac
    dat['flux_extcorr'] = f_extcorr
    dat['unc_extcorr'] = ef_extcorr

    dat.loc[isdet, 'mag'] = -2.5*np.log10(fdet*1E-6)+8.90
    dat.loc[isdet, 'emag'] = (2.5/np.log(10)) * (efdet/fdet)
    fdet_corr = dat['flux_extcorr'][isdet]
    dat['mag_extcorr'] = np.copy(dat['mag'])
    dat.loc[isdet, 'mag_extcorr'] = -2.5*np.log10(fdet_corr*1E-6)+8.90

    # Non-detections / upper limits
    # Calculate an upper limit (limiting magnitude) for all exposures
    dat['maglim'] = -2.5*np.log10(dat['dujy']*1E-6*SNU)+8.90
    dat['maglim_extcorr'] = -2.5*np.log10(dat['unc_extcorr']*1E-6*SNU)+8.90

    # Now create a new dataframe with the same headings as our other dfs
    d = {}
    d['#instrument'] = np.array(['PS1 '+v for v in dat['pscamera'].values])
    d['mjdstart'] = dat['#mjd'].values
    d['exp'] = dat['exptime'].values
    d['flt'] = dat['filter'].values
    d['flux'] = dat['ujy'].values
    d['flux_extcorr'] = dat['flux_extcorr'].values
    d['unc'] = dat['dujy'].values
    d['unc_extcorr'] = dat['unc_extcorr'].values
    d['sig'] = np.abs(dat['ujy'].values/dat['dujy'].values)
    d['flare'] = np.logical_and(d['sig']>5, d['mjdstart']>59856.4)
    d['mag'] = dat['mag'].values
    d['mag_extcorr'] = dat['mag_extcorr'].values
    d['emag'] = dat['emag'].values
    d['maglim'] = dat['maglim'].values
    d['maglim_extcorr'] = dat['maglim_extcorr'].values

    df = pd.DataFrame(d)

    # Keep only the dates after 2022-08-25
    keep = df['mjdstart'] > Time("2022-08-25T00:00:00", format='isot').mjd
    df = df[keep]

    return df



def get_chimera():
    """ Get CHIMERA data """
    inputf = ddir + "/lc_chimera.txt"
    dat  = pd.read_fwf(inputf)

    # Add a magnitudes column
    dat['mag'] = [99]*len(dat)
    dat['emag'] = [99]*len(dat)
    dat['maglim'] = [99]*len(dat)

    # Detections
    SNU = 3 # provide 3-sigma U.L.
    SNT = 3
    isdet = dat['sig']>=SNT # confident detection
    fdet = dat['flux'][isdet]
    efdet = dat['unc'][isdet]

    # Correct for MW extinction
    f_extcorr = np.copy(dat['flux'].values)
    ef_extcorr = np.copy(dat['unc'].values)
    for f in np.unique(dat['flt']):
        choose = dat['flt'].values==f
        fac = 10**(vals.ext[f]/2.5)
        f_extcorr[choose] = dat['flux'].values[choose]*fac
        ef_extcorr[choose] = dat['unc'].values[choose]*fac
    dat['flux_extcorr'] = f_extcorr
    dat['unc_extcorr'] = ef_extcorr

    dat.loc[isdet, 'mag'] = -2.5*np.log10(fdet*1E-6)+8.90
    dat.loc[isdet, 'emag'] = (2.5/np.log(10)) * (efdet/fdet)
    fdet_corr = dat['flux_extcorr'][isdet]
    dat['mag_extcorr'] = np.copy(dat['mag'])
    dat.loc[isdet, 'mag_extcorr'] = -2.5*np.log10(fdet_corr*1E-6)+8.90

    # Non-detections / upper limits
    # Calculate an upper limit (limiting magnitude) for all exposures
    dat['maglim'] = -2.5*np.log10(dat['unc']*1E-6*SNU)+8.90
    dat['maglim_extcorr'] = -2.5*np.log10(dat['unc_extcorr']*1E-6*SNU)+8.90

    return dat



def get_full_opt():
    """ Retrieve a table of ALL optical photometry 
    Report upper limits as 3-sigma
    Label flares as 5-sigma
    MJD is the start time
    """
    # Put the ZTF data into a dictionary
    jd,exp,filt,mag,mag_extcorr,emag,fujy,efujy,fujy_extcorr,efujy_extcorr = \
            get_ipac() # ULs are 3-sig

    add_dict = {}
    add_dict['#instrument'] = ['ZTF']*len(jd)
    add_dict['mjdstart'] = Time(jd, format='jd').mjd
    add_dict['exp'] = exp
    add_dict['flt'] = filt
    add_dict['flux'] = fujy
    add_dict['unc'] = efujy
    add_dict['flux_extcorr'] = fujy_extcorr
    add_dict['unc_extcorr'] = efujy_extcorr
    add_dict['sig'] = np.abs(fujy/efujy)
    add_dict['mag'] = mag
    add_dict['mag_extcorr'] = mag_extcorr
    add_dict['emag'] = emag

    # Give all observations a 3-sigma limiting magnitude
    add_dict['maglim'] = -2.5*np.log10(efujy*1E-6*3)+8.90
    add_dict['maglim_extcorr'] = -2.5*np.log10(efujy_extcorr*1E-6*3)+8.90
    add_dict = pd.DataFrame(add_dict)

    # Add Dan's photometry to the ZTF photometry
    dat = get_dan_lc() # 3-sigma
    ztf_dan = dat.append(add_dict, ignore_index=True)

    # Add Dan's CHIMERA photometry
    chimera = get_chimera()
    ztf_dan_chimera = ztf_dan.append(chimera, ignore_index=True)

    # Add the PS1 photometry
    ps1 = get_panstarrs()
    ztf_dan_chimera_ps1 = ztf_dan_chimera.append(ps1, ignore_index=True)

    # Add the Lulin photometry
    lulin = get_lulin() # 3-sigma
    lulin_tel = []
    for f in lulin['filt'].values:
        if f=='g':
            lulin_tel.append('LOT')
        elif f=='r':
            lulin_tel.append('SLT')
        else:
            print("problem")
    add_dict = {}
    add_dict['#instrument'] = lulin_tel
    add_dict['mjdstart'] = lulin['MJD']
    add_dict['exp'] = ['300']*len(lulin) # asking group
    add_dict['flt'] = lulin['filt']
    add_dict['flux'] = '' # not provided
    add_dict['unc'] = '' # not provided
    add_dict['sig'] = [-99]*len(lulin) # not provided
    add_dict['flare'] = ['']*len(lulin) # none
    add_dict['mag'] = [99]*len(lulin) # not provided
    add_dict['emag'] = [99]*len(lulin) # not provided
    add_dict['maglim'] = lulin['maglim']

    # Correct for MW extinction
    maglim_extcorr = []
    for i,b in enumerate(lulin['filt'].values):
        maglim_extcorr.append(lulin['maglim'].values[i]-vals.ext[b])
    add_dict['maglim_extcorr'] = maglim_extcorr

    add_dict = pd.DataFrame(add_dict)

    full_dict = ztf_dan_chimera_ps1.append(add_dict, ignore_index=True)

    # Indicate that all rows of this table are single images
    full_dict['nobs'] = [1]*len(full_dict)

    # Mark flares
    full_dict['isflare'] = np.logical_and(full_dict['mjdstart']>59856.4,
                                          full_dict['sig']>5)

    # Mark the transient
    full_dict['istransient'] = np.logical_and(full_dict['mjdstart']<59856.4,
                                          full_dict['sig']>3)

    return full_dict.sort_values(by=['mjdstart'], ignore_index=True, axis=0)


def get_dan_lc():
    """ Get the LC data provided by Dan Perley 
    Upper limits are 3-sigma """
    inputf = ddir + "/full_lc.txt"
    dat  = pd.read_fwf(inputf, infer_nrows=200)

    # Add a magnitudes column
    dat['mag'] = [99]*len(dat)
    dat['emag'] = [99]*len(dat)
    dat['maglim'] = [99]*len(dat)

    # Detections
    SNU = 3 # provide 3-sigma U.L.
    SNT = 3
    isdet = dat['sig']>=SNT # confident detection
    fdet = dat['flux'][isdet]
    efdet = dat['unc'][isdet]

    # Correct for MW extinction
    f_extcorr = np.copy(dat['flux'].values)
    ef_extcorr = np.copy(dat['unc'].values)
    for f in np.unique(dat['flt']):
        choose = dat['flt'].values==f
        fac = 10**(vals.ext[f]/2.5)
        f_extcorr[choose] = dat['flux'].values[choose]*fac
        ef_extcorr[choose] = dat['unc'].values[choose]*fac
    dat['flux_extcorr'] = f_extcorr
    dat['unc_extcorr'] = ef_extcorr

    dat.loc[isdet, 'mag'] = -2.5*np.log10(fdet*1E-6)+8.90
    dat.loc[isdet, 'emag'] = (2.5/np.log(10)) * (efdet/fdet)
    fdet_corr = dat['flux_extcorr'][isdet]
    dat['mag_extcorr'] = np.copy(dat['mag'])
    dat.loc[isdet, 'mag_extcorr'] = -2.5*np.log10(fdet_corr*1E-6)+8.90

    # Non-detections / upper limits
    # Calculate an upper limit (limiting magnitude) for all exposures
    dat['maglim'] = -2.5*np.log10(dat['unc']*1E-6*SNU)+8.90
    dat['maglim_extcorr'] = -2.5*np.log10(dat['unc_extcorr']*1E-6*SNU)+8.90

    return dat

