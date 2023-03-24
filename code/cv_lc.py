""" Check to see if any dwarf nova light curves resemble AT2022tsd """

import extinction
from dustmaps.sfd import SFDQuery
from dustmaps.config import config
config['data_dir'] = '/Users/annaho/Dropbox/astro/tools/scanning/sfddata-master/'
from penquins import Kowalski
from astropy.time import Time
from astropy.coordinates import SkyCoord
import random
import vals
import itertools
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from helpers import *


def read_ipac_forced_phot(inputf):
    """ Get the full IPAC forced phot light curve """
    a = pd.read_table(inputf, comment='#', delimiter=' ')

    # Extract all field/ccd/filter combinations
    fid = np.zeros(len(a['filter,']))
    fid[a['filter,']=='ZTF_g'] = 1
    fid[a['filter,']=='ZTF_r'] = 2
    fid[a['filter,']=='ZTF_i'] = 3
    fcqf = (a['field,']*10000+a['ccdid,']*100+a['qid,']*10+fid).astype(int).values
    ufcqf = np.unique(fcqf)

    # Bad pixel contamination
    # leaving this out for now...
    if a['procstatus'].dtype=='int64':
        bad_pix = a['procstatus'].values == 56
    else:
        bad_pix = np.logical_or(a['procstatus'].values == '56',
            a['procstatus'].values == '56,57')

    # Bad observing conditions
    bad_conditions = np.logical_or.reduce((
            a['scisigpix,'].values > 5*np.median(a['scisigpix,'].values),
            a['zpmaginpsci,'].values > 5*np.median(a['zpmaginpsci,'].values),
            a['zpmaginpscirms,'].values > 5*np.median(a['zpmaginpscirms,'].values)))

    # Get the full light curve
    jd = a['jd,'].values
    flux = a['forcediffimflux,'].values
    eflux = a['forcediffimfluxunc,'].values
    filt = a['filter,'].values
    maglim = a['diffmaglim,'].values
    zp = a['zpdiff,'].values
    programid = a['programid,'].values
    refstart = a['refjdstart,'].values
    refend = a['refjdend,'].values

    # some flux values are NaN
    nan = np.isnan(flux)

    # some chisq values are bad
    bad_chisq = np.isnan(a['forcediffimchisq,'])

    # Discard bad values
    #good = np.logical_and.reduce((~bad_pix, ~bad_conditions, ~nan, ~bad_chisq))
    good = np.logical_and.reduce((~bad_conditions, ~nan, ~bad_chisq))
    jd = jd[good]
    flux = flux[good]
    eflux = eflux[good]
    filt = filt[good]
    maglim = maglim[good]
    zp = zp[good]
    fcqf = fcqf[good]
    programid = programid[good]
    refstart = refstart[good]
    refend = refend[good]
    ufcqf = np.unique(fcqf)

    # Cobble together the LC
    SNT = 3
    SNU = 5

    order = np.argsort(jd)
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
    filt = filt_final
    refstart = refstart[order]
    refend = refend[order]

    return jd,filt,fujy,efujy,mag,emag,refstart,refend


def get_forced_phot_lc():
    """ 
    Get forced photometry LC for the 182 CVs from BTS with mpeak > 18 mag """
    cvs = pd.read_csv("../data/btscvs.txt")
    names = list(cvs['ZTFID'].values)

    # Log onto Kowalski
    s = logon()
    
    # Get IPAC command for each object
    print("Name", "RA", "Dec", "ForcePhot")
    for name in names:
        #print(name)
        ra,dec = get_pos(s, name)
        print(name,ra,dec)
        # JD start and end correspond to March 2018 and 1 Jan 2023
        cmd = 'wget --http-user=ztffps --http-passwd=dontgocrazy! -O log.txt "https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?ra=%s&dec=%s&jdstart=2458178&jdend=2459945.5&email=annayqho@gmail.com&userpass=hasa654"' %(ra,dec)
        print(cmd)


def analyze_lc():
    """ 
    Analyze the force phot LC to look for short-timescale variability """
    cvs = pd.read_csv("cv_table.txt", delimiter=' ')
    
    for j,cv in enumerate(cvs['Name'].values):
        print(j)
        # Construct filename
        f = '../data/opt/ipac_cv/%s.txt' %cv
        jd,filt,fujy,efujy,mag,emag,refstart,refend = \
                read_ipac_forced_phot(f)
        # Grab obs that have another obs w/in the same night
        jd_int = jd.astype(int)
        jd_int_unique, counts = np.unique(jd_int, return_counts=True)
        keep_jd = counts > 1
        keep = np.array([jdval in jd_int_unique[keep_jd] for jdval in jd_int])
        jd = jd[keep]
        jd_int = jd.astype(int)
        filt = filt[keep]
        fujy = fujy[keep]
        efujy = efujy[keep]
        mag = mag[keep]
        emag = emag[keep]
        refstart = refstart[keep]
        refend = refend[keep]
        # For each JD, calculate the baseline
        names = []
        jdvals = []
        baseline = []
        dflux = []
        edflux = []
        for jdval in np.unique(jd_int):
            use = jd_int==jdval
            # Get the observations for that JD 
            jd_use = jd[use]
            filt_use = filt[use]
            fujy_use = fujy[use]
            efujy_use = efujy[use]
            mag_use = mag[use]
            emag_use = emag[use]
            refstart_use = refstart[use]
            refend_use = refend[use]
            # Check that the minimum difference is < 0.7d
            diffs = jd_use[1:]-jd_use[0:-1]
            if min(diffs)<0.7: # then continue
                # Check that there is at least one detection
                sig = fujy_use / efujy_use
                has_det = sum(sig>5)>0
                if has_det: # then continue
                    # Check that there is a significant change in flux over < 0.7d
                    # Get all pairs of JD, Flux, eFlux
                    jd_pairs = list(itertools.combinations(jd_use, 2))
                    f_pairs = list(itertools.combinations(fujy_use, 2))
                    ef_pairs = list(itertools.combinations(efujy_use, 2))
                    filt_pairs = list(itertools.combinations(filt_use, 2))
                    refstart_pairs = list(itertools.combinations(refstart_use, 2))
                    refend_pairs = list(itertools.combinations(refend_use, 2))
                    for i,jd_pair in enumerate(jd_pairs):
                        # only bother with the same filter
                        filt_pair = filt_pairs[i]
                        same_filt = filt_pair[0]==filt_pair[1]
                        # require same ref start (to get rid of really variable things)
                        refstart_pair = refstart_pairs[i]
                        refend_pair = refend_pairs[i]
                        same_ref = np.logical_and(
                            refstart_pair[0]==refstart_pair[1],
                            refend_pair[0]==refend_pair[1])
                        dt = jd_pair[1]-jd_pair[0]
                        if np.logical_and.reduce((same_filt, same_ref, dt < 0.7)): 
                            jdvals.append(jdval)
                            baseline.append(dt)
                            f_pair = f_pairs[i]
                            ef_pair = ef_pairs[i]

                            # Look for points with a change of > an OOM
                            fratio = np.abs(f_pair[1]/f_pair[0])
                            high_amp = np.logical_or(fratio>10, fratio<0.1)
                            # and have significantly different fluxes
                            df = np.abs(f_pair[1]-f_pair[0])
                            edf = np.sqrt(ef_pair[1]**2+ef_pair[0]**2)
                            sig = df/edf > 3
                            # combined requirements
                            if np.logical_and(high_amp, sig):
                                print(cv, jdval)
                                print(fratio)
                                print(f_pair)
                                print(filt_pair)
                                t0 = min(jd_use)
                                fs = ['g', 'r', 'i']
                                col = ['Aquamarine', 'Crimson', 'Goldenrod']
                                for k,f in enumerate(fs):
                                    choose = filt_use==f
                                    if sum(choose)>0:
                                        plt.errorbar(
                                            (jd_use[choose]-t0)*24*60, 
                                            fujy_use[choose], 
                                            efujy_use[choose], 
                                            fmt='o', c=col[k])
                                plt.xlabel("Minutes")
                                plt.ylabel("Flux")
                                #plt.show()
                                plt.title("filt pair: %s, %s" %(filt_pair[0], filt_pair[1]))
                                plt.tight_layout()
                                plt.savefig("sig_dflux_%s_%s.png" %(cv,jdval))
                                plt.close()
                                dflux.append(df)
                                edflux.append(edf)
                                break
                                

        #baseline = np.array(baseline)
        #dflux = np.array(dflux)
        #edflux = np.array(edflux)
        #if max(np.abs(dflux/edflux))>3:
        #    print(cv)


if __name__=="__main__":
    analyze_lc()
