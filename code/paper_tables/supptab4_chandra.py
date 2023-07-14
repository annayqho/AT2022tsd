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
            ['$t_\mathrm{start}$', '$\Delta t$\\footnote{Rest frame}', 
             '$t_\mathrm{exp}$', '$F_X$', '$L_X$'])
    unit_headings = np.array(
            ['(UT)', '(days)', '(ksec)', 
             '($10^{-14}$ erg\,s$^{-1}$\,cm$^{-2}$)', 
             '$(10^{43}$\,erg\,s$^{-1}$)'])
    label = "tab:chandra"

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

    caption="Chandra X-ray Observatory observations of AT2022tsd."

    outputf = open("paper_table_%s.txt" %label, "w")
    outputf.write("\\begin{table*} \n")
    outputf.write("\\centering \n")
    outputf.write("\\begin{tabular}{%s}\n" %colstr)
    outputf.write("\hline\hline\n")
    outputf.write(colheadstr+'\\\ \n')
    outputf.write(unitstr+'\\\ \n')
    outputf.write("\hline\n")

    # Get data
    df = load_chandra()

    for i in np.arange(len(df)):
        # Convert date to readable date and time
        t = Time(df['MJD'][i], format='mjd')
        tstr = t.isot.replace("T", " ")[0:16]

        # Rest-frame days
        dtstr = '{:.2f}'.format((t.jd-vals.t0)/(1+vals.z))

        # Exposure time
        texp = int(df['Exp'][i])

        # Flux
        f = '{:.2f}'.format(df['Flux'][i]/1E-14)
        # Errors are asymmetric
        ef_up = '{:.2f}'.format(df['uFlux'][i]/1E-14)
        ef_bot = '{:.2f}'.format(df['lFlux'][i]/1E-14)
        fstr = "$%s^{+%s}_{-%s}$" %(f,ef_up,ef_bot)

        # Luminosity
        L = '{:.2f}'.format(df['L'][i]/1E43)
        eL_up = '{:.2f}'.format(df['Lpos'][i]/1E43)
        eL_bot = '{:.2f}'.format(df['Lneg'][i]/1E43)
        Lstr = "$%s^{+%s}_{-%s}$" %(L,eL_up,eL_bot)

        row = rowstr %(tstr,dtstr,texp,fstr,Lstr)
        print(row)
        outputf.write(row)

    outputf.write("\hline \n")
    outputf.write("\\end{tabular} \n")
    outputf.write("\\caption{%s} \n" %caption)
    outputf.write("\\label{%s} \n" %label)

    outputf.write("\end{table*} \n")
    outputf.close()


if __name__=="__main__":
    print_table()
