""" Print a table of optical photometry """
import numpy as np
import sys
from astropy.time import Time
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
from get_opt import get_ipac
from get_radio_at2022tsd import *
from get_xray import *
import vals

def print_table():
    """ Print the table """

    # Headings
    headings = np.array(
            ['t', 't_\mathrm{err}', '$\Delta t$\\footnote{Rest frame}', 
             '$\Delta T$', 'Count Rate', '$F_X$', '$L_X$', 'Telescope'])
    unit_headings = np.array(
            ['(UT)', '(days)', '(days)', '(ksec)', '($10^{-3}$\,s$^{-1}$)', 
             '$10^{-14}$ (erg\,s$^{-1}$\,cm$^{-2}$)', '$(10^{43}$\,erg\,s$^{-1}$)', ''])
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

    # Error bars are asymmetric
    df = load_swift()

    for i in np.arange(len(t)):
        # Convert date to readable date and time
        t = Time(df['!MJD    '][i], format='mjd')
        et = df['T_+ve   '][i]

        # Effective time...only want minute precision
        tstr = str((t-et).isot).replace('T', ' ').split('.')[0][0:-3]
        # Error on the time
        terror = '{:.2f}'.format(et)
        # Rest-frame days
        dtstr = '{:.2f}'.format((t.jd-vals.t0)/(1+vals.z))
        # Exposure time
        texp = '{:.2f}'.format(get_exp(i+1)/1000)

        # Counts
        ct = '{:.2f}'.format(df['Rate    '][i]/1E-3)
        # Errors are symmetric
        ect = '{:.2f}'.format(df['Ratepos '][i]/1E-3)
        ctstr = "$%s\pm%s$" %(ct,ect)

        # Flux
        f = '{:.2f}'.format(df['Flux'][i]/1E-14)
        # Errors are symmetric
        ef = '{:.2f}'.format(df['Fluxpos'][i]/1E-14)
        fstr = "$%s\pm%s$" %(ct,ect)

        # Luminosity
        L = '{:.2f}'.format(df['L'][i]/1E43)
        # Errors are symmetric
        eL = '{:.2f}'.format(df['Lpos'][i]/1E43)
        Lstr = "$%s\pm%s$" %(ct,ect)

        row = rowstr %(tstr,terror,dtstr,texp,ctstr,fstr,Lstr,'Swift/XRT')
        print(row)
        outputf.write(row)

    outputf.write("\hline \n")

    outputf.write("\end{longtable} \n")
    outputf.write("\end{center} \n")
    outputf.close()


if __name__=="__main__":
    print_table()
