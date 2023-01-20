""" Print a table of optical photometry """
import numpy as np
import sys
from astropy.time import Time
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
from get_opt import get_ipac
import values

def print_table():
    """ Print the table """

    # Headings
    headings = np.array(
            ['Date', '$\Delta t$\\footnote{Rest frame}', 'Filter', 
             'Mag\\footnote{Not corrected for Galactic extinction. Upper limits are 5-$\sigma$.}', 
             'eMag', 'Instrument'])
    unit_headings = np.array(
            ['(UT)', '(days)', '', 
             '(AB)', '(AB)', ''])
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

    jd,filt,mag,emag = get_ipac()

    for lc_point in list(zip(jd,filt,mag,emag)):
        # Convert JD to readable date and time
        tstr = Time(lc_point[0], format='jd').isot.replace('T', ' ').split('.')[0]
        dtstr = '{:.4f}'.format((lc_point[0]-values.t0)/(1+values.z))
        filtstr = "$\mathrm{ZTF}_{%s}$" %lc_point[1]
        # Upper limit
        if lc_point[3]==99:
            mstr = '$>{:.2f}$'.format(lc_point[2])
            emstr = '--'
        # Detection
        else:
            mstr = '${:.2f}$'.format(lc_point[2])
            emstr = '${:.2f}$'.format(lc_point[3])
        row = rowstr %(tstr,dtstr,filtstr,mstr,emstr,'ZTF')
        print(row)
        outputf.write(row)

    outputf.write("\hline \n")

    outputf.write("\end{longtable} \n")
    outputf.write("\end{center} \n")
    outputf.close()


if __name__=="__main__":
    print_table()
