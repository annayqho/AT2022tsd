""" Print a table of optical photometry """
import numpy as np
import sys
from astropy.time import Time
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
from get_opt import get_ipac
from get_radio_at2022tsd import *
from get_xray import load_swift
import vals

def print_table():
    """ Print the table """

    # Headings
    headings = np.array(
            ['Start Date', '$\Delta t$\\footnote{Rest frame}', 'Count Rate', 
             '$F_X$', '$L_X$', 'Telescope'])
    unit_headings = np.array(
            ['(UT)', '(days)', '(s$^{-1}$)', 
             '(erg\,s$^{-1}$\,cm$^{-2}$)', '(erg\,s$^{-1}$)', ''])
    label = "tab:xray-observations"

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

    caption="X-ray observations of AT2022tsd."

    outputf = open("paper_table_%s.txt" %label, "w")
    outputf.write("\\begin{center} \n")
    outputf.write("\\begin{longtable}{%s} \n" %colstr)
    outputf.write("\\caption{%s} \n" %caption)
    outputf.write("\\label{%s}\\\ \n" %label)
    outputf.write("\hline\hline\n")
    outputf.write(colheadstr+'\\\ \n')
    outputf.write(unitstr+'\\\ \n')
    outputf.write("\hline\n")

    t,Ls,lLs,uLs = load_swift()

    for i in np.arange(len(dat['Date'])):
        # Convert date to readable date and time
        tstr = dat['Date'].values[i].replace('T', ' ').split('.')[0]
        # Rest-frame days
        dtstr = '{:.2f}'.format((Time(dat['Date'].values[i]).jd-vals.t0)/(1+vals.z))
        nustr = int(dat['Freq_Obs'].values[i])
        # I think they're all detections
        fstr = '${:.3f}$'.format(dat['Flux'].values[i])
        efstr = '${:.3f}$'.format(dat['eFlux'].values[i])
        row = rowstr %(tstr,dtstr,nustr,fstr,efstr,dat['Tel'].values[i])
        print(row)
        outputf.write(row)

    outputf.write("\hline \n")

    outputf.write("\end{longtable} \n")
    outputf.write("\end{center} \n")
    outputf.close()


if __name__=="__main__":
    print_table()
