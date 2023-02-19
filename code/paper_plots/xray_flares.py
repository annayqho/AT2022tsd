""" Plot the Chandra X-ray flares """

import pandas as pd
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt

dd = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/xray"
fs = ['26641/repro/src2_sub_lc_500s.fits','26642/repro/src2_sub_lc_500s.fits',
      '26643/repro/src2_sub_lc_500s.fits','26644/repro/src2_sub_lc_500s.fits',
      '26645/repro/src2_sub_lc_500s.fits','27639/repro/src2_sub_lc_500s.fits',
      '27643/repro/src2_sub_lc_500s.fits']

f = fs[0]
tab = pyfits.open(dd+"/"+f)
dt = tab.get_column('dt').values
rate = tab.get_column('net_rate').values
erate = tab.get_column('err_rate').values
plt.errorbar(
        dt, rate, yerr=erate, marker="o", color="red", mfc="black",
        mec="black", ecolor="grey")


plt.show()

