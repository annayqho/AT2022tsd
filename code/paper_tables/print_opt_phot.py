""" Print a table of optical photometry """
import numpy as np
import pandas as pd
import sys
from astropy.time import Time
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
from get_opt import get_ipac,get_flares
import vals

def print_table():
    """ Print the table """

    # Headings
    headings = np.array(
            ['Date', '$\Delta t$\\footnote{Rest frame}', 'Filter', 
             'Mag\\footnote{Not corrected for Galactic extinction. Upper limits are 5-$\sigma$.}', 
             'eMag', 'Instrument', 'Flare?\\footnote{$>5$-$\sigma$ detections.}'])
    unit_headings = np.array(
            ['(UT)', '(days)', '', 
             '(AB)', '(AB)', '', ''])
    label = "optical-photometry"

    ncol = len(headings)
    colstr = ""
    colstr += 'l'
    for col in np.arange(ncol-1): colstr+="l"
    print(colstr)

    colheadstr = ""
    for col in np.arange(ncol-1):
        colheadstr += "%s & " %headings[col]
    colheadstr += "%s" %headings[-1]

    unitstr = ""
    for col in np.arange(ncol-1):
        unitstr += "%s & " %unit_headings[col]
    unitstr += "%s" %unit_headings[-1]

    rowstr = ""
    for col in np.arange(ncol-1):
        rowstr += "%s & "
    rowstr += "%s \\\ \n"

    caption="Optical photometry of AT2022tsd."

    outputf = open("paper_table_%s.txt" %label, "w")
    outputf.write("\\begin{center} \n")
    outputf.write("\\begin{longtable}{%s} \n" %colstr)
    outputf.write("\\caption{%s} \n" %caption)
    outputf.write("\\label{%s}\\\ \n" %label)
    outputf.write("\hline\hline\n")
    outputf.write(colheadstr+'\\\ \n')
    outputf.write(unitstr+'\\\ \n')
    outputf.write("\hline\n")

    # The ZTF lines
    jd,exp,filt,mag,emag,fujy,efujy = get_ipac()
    tel = np.array(['ZTF']*len(jd))
    flare = np.array(['']*len(jd))
    limmag = np.copy(mag)

    # the IPAC flares are from JD > 2459857
    is_flare = np.logical_and(emag<99, jd>2459854)
    flare[is_flare] = ['*']*sum(is_flare)

    # The flare lines
    tel_fl,mjd_fl,filt_fl,mag_fl,emag_fl,limmag_fl,flare_fl = get_flares()
    flare_fl[pd.isnull(flare_fl)] = np.array(['']*sum(pd.isnull(flare_fl)))

    # Combine
    jd = np.hstack((jd, Time(mjd_fl, format='mjd').jd))
    filt = np.hstack((filt, filt_fl))
    mag = np.hstack((mag, mag_fl))
    emag = np.hstack((emag, emag_fl))
    flare = np.hstack((flare, flare_fl))
    tel = np.hstack((tel, tel_fl))
    limmag = np.hstack((limmag, limmag_fl))

    # Sort
    order = np.argsort(jd)
    jd = jd[order]
    filt = filt[order]
    mag = mag[order]
    emag = emag[order]
    flare = flare[order]
    tel = tel[order]
    limmag = limmag[order]

    for i in np.arange(len(jd)):
        # Convert JD to readable date and time
        tstr = Time(jd[i], format='jd').isot.replace('T', ' ').split('.')[0]

        # Calculate the dt in the rest frame
        dtstr = '{:.4f}'.format((jd[i]-vals.t0)/(1+vals.z))

        # Filter
        filtstr = filt[i]
        if tel[i]=='ZTF':
            filtstr = "$\mathrm{ZTF}_{%s}$" %filt[i]

        # Upper limit
        is_nondet = np.logical_or.reduce(
                (emag[i]==99, np.isnan(emag[i]), mag[i]>limmag[i]))
        if is_nondet:
            mstr = '$>{:.2f}$'.format(limmag[i])
            emstr = '--'
        # Detection
        else:
            mstr = '${:.2f}$'.format(mag[i])
            emstr = '${:.2f}$'.format(emag[i])

        # Flare string
        flstr = flare[i]
        if flstr=='nan':
            flstr = ''

        row = rowstr %(tstr,dtstr,filtstr,mstr,emstr,tel[i],flstr)
        print(row)
        outputf.write(row)

    outputf.write("\hline \n")

    outputf.write("\end{longtable} \n")
    outputf.write("\end{center} \n")
    outputf.close()


if __name__=="__main__":
    print_table()
