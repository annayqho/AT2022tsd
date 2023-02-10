""" Code to get the optical photometry """

import pandas as pd
import numpy as np
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
    SNU = 5 # S/N to use when computing an upper limit

    mag = np.array([99]*len(flux)).astype(float)
    emag = np.array([99]*len(eflux)).astype(float)
    is_det = flux/eflux > SNT
    mag[is_det] = zp[is_det]-2.5*np.log10(flux[is_det])
    emag[is_det] = 1.0857*eflux[is_det]/flux[is_det]
    mag[~is_det] = zp[~is_det]-2.5*np.log10(SNU*eflux[~is_det])

    # Convert the flux and eflux into physical units, uJy
    f0 = 10**(0.4*zp)
    fratio = flux/f0
    fujy = fratio * 3631 * 1E6
    efujy = eflux * (fujy/flux)

    order = np.argsort(jd)
    filt_final = np.copy(filt[order])
    filt_final[filt_final=='ZTF_g'] = 'g'
    filt_final[filt_final=='ZTF_r'] = 'r'
    filt_final[filt_final=='ZTF_i'] = 'i'

    # Final arrays for the light curve
    jd = jd[order]
    fujy = fujy[order]
    efujy = efujy[order]
    mag = mag[order]
    emag = emag[order]
    exp = exp[order]
    filt = np.copy(filt_final)

    return jd,exp,filt,mag,emag,fujy,efujy


def get_transient():
    """ Get all the opt phot of the main transient event """
    jd,filt,mag,emag,fujy,efujy = get_ipac()


def get_ultraspec():
    """ Get the ULTRASPEC observations """
    inputf = ddir + "/flares_lris_ultraspec.txt"
    dat  = pd.read_fwf(
            inputf, comment='#', 
            names=['MJD','Exp','Filt','Flux','Unc'])
    return dat


def get_non_ztf():
    # This is everything except ZTF.
    inputf = ddir + "/flares_lris_ultraspec_lt.txt"
    dat  = pd.read_fwf(inputf, comment='#', 
                         names=['Tel','MJD','Exp','Filt','Flux','Unc'])
    # Convert to mag
    dat['Mag'] = [99]*len(dat)
    dat['eMag'] = [99]*len(dat)
    dat['Maglim'] = [99]*len(dat)
    dat['Flare'] = ['']*len(dat)
    SNU = 5
    SNT = 3
    SNFlare = 5
    isdet = dat['Flux']/dat['Unc']>SNT # confident detection
    fdet = dat['Flux'][isdet]
    efdet = dat['Unc'][isdet]
    dat.loc[isdet, 'Mag'] = -2.5*np.log10(fdet*1E-6)+8.90
    dat.loc[isdet, 'eMag'] = (2.5/np.log(10)) * (efdet/fdet)
    dat.loc[~isdet, 'Maglim'] = -2.5*np.log10(
            dat['Unc'][~isdet]*1E-6*SNU)+8.90
    dat.loc[dat['Flux']/dat['Unc']>SNFlare, 'Flare'] = ['*']*sum(
            dat['Flux']/dat['Unc']>SNFlare)
    return dat
