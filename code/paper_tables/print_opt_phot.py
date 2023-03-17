""" Print a table of optical photometry """
import numpy as np
import pandas as pd
import sys
from astropy.time import Time
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
from get_opt import get_full_opt,get_last
import vals

def print_table_all():
    """ Print the table """

    # Headings
    headings = np.array(
            ['Date', '$\Delta t$\\footnote{Rest frame}', '$t_\mathrm{exp}$', 
             'Filter', 
             'Mag\\footnote{Not corrected for Galactic extinction.}', 
             'eMag\\footnote{Upper limits reported as 3-$\sigma$.}', 
             'Instrument', 'Flare?\\footnote{$>5$-$\sigma$ detections.}'])
    unit_headings = np.array(
            ['(UT)', '(days)', '(sec.)', '',
             '(AB)', '(AB)', '', ''])
    label = "tab:optical-photometry"

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

    # Get the data
    dat = get_full_opt()

    for i in np.arange(len(dat)):
        # Get the things you need
        mjd = dat['mjdstart'][i]
        jd = Time(mjd, format='mjd').jd
        tel = dat['#instrument'][i]
        filt = dat['flt'][i]
        exptime = dat['exp'].values[i]
        mag = dat['mag'][i]
        emag = dat['emag'][i]
        limmag = dat['maglim'][i]
        isflare = dat['isflare'][i]

        # Convert to readable date and time
        tstr = Time(
                mjd, format='mjd').isot.replace('T', ' ').split('.')[0]

        # Calculate the dt in the observer frame
        dtstr = '{:.4f}'.format(jd-vals.t0)

        # Filter
        filtstr = filt
        if tel=='ZTF':
            filtstr = "${%s}_\mathrm{ZTF}$" %filt
        else:
            filtstr = "$%s$" %filt

        # Upper limit
        is_nondet = np.logical_or.reduce(
                (emag==99, np.isnan(emag), mag>limmag))
        if is_nondet:
            mstr = '$>{:.2f}$'.format(limmag)
            emstr = '--'
        # Detection
        else:
            mstr = '${:.2f}$'.format(mag)
            emstr = '${:.2f}$'.format(emag)

        # Flare string
        flstr = ''
        if isflare:
            flstr = '*'

        row = rowstr %(tstr,dtstr,exptime,filtstr,mstr,emstr,tel,flstr)
        print(row)
        outputf.write(row)

    outputf.write("\hline \n")

    outputf.write("\end{longtable} \n")
    outputf.write("\end{center} \n")
    outputf.close()


def print_table_last():
    """ Print the table """
    # Headings
    headings = np.array(
            ['Date', '$\Delta t$\\footnote{Observer frame}', 
             'Flux\\footnote{Not corrected for Galactic extinction.\
                     ZP=25 (AB), calibrated to Gaia $G_p$ band.}',\
             'eFlux', '\# Frames'])
    unit_headings = np.array(['(UT)', '(days)', '', '', ''])
    label = "tab:last-photometry"

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

    caption="LAST observations of AT2022tsd, binned into two-minute exposures."

    outputf = open("paper_table_%s.txt" %label, "w")
    outputf.write("\\begin{center} \n")
    outputf.write("\\begin{longtable}{%s} \n" %colstr)
    outputf.write("\\caption{%s} \n" %caption)
    outputf.write("\\label{%s}\\\ \n" %label)
    outputf.write("\hline\hline\n")
    outputf.write(colheadstr+'\\\ \n')
    outputf.write(unitstr+'\\\ \n')
    outputf.write("\hline\n")

    # Get data
    mjd,flux,eflux,nobs = get_last()

    for i in np.arange(len(mjd)):
        # Convert to readable date and time
        tstr = Time(
                mjd[i], format='mjd').isot.replace('T', ' ').split('.')[0]

        # Calculate the dt in the rest frame
        jd = Time(mjd[i], format='mjd').jd
        dtstr = '{:.4f}'.format((jd-vals.t0)/(1+vals.z))

        fstr = '{:.2f}'.format(flux[i])
        efstr = '{:.2f}'.format(eflux[i])
        nobsstr = nobs[i]

        row = rowstr %(tstr,dtstr,fstr,efstr,nobsstr)
        print(row)
        outputf.write(row)

    outputf.write("\hline \n")

    outputf.write("\end{longtable} \n")
    outputf.write("\end{center} \n")
    outputf.close()


if __name__=="__main__":
    print_table_all() # the main table
    #print_table_last() # the table of LAST photometry
