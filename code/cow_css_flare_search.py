""" Analyze the CSS161010 forced photometry from ATLAS """

import pandas as pd

dat = pd.read_fwf("css161010_atlas_forced_phot.txt")
t = dat['###MJD'].values
dt = t-t[0]
f = dat['uJy'].values
ef = dat['duJy'].values
limmag_raw = dat['mag5sig Sky'].values
limmag = np.array([val.split(' ')[0] for val in limmag_raw]).astype(float)

choose = np.logical_and(dt>20, dt<100)
dt = dt[choose]
f = f[choose]
ef= ef[choose]
limmag = limmag[choose]

choose = limmag >= 19.2
dt = dt[choose]
dt_hr = dt*24
dt_hr_int = np.unique((dt_hr.astype(int)))

# How many images have dt < 100 and limmag < 19.2
choose = np.logical_and.reduce((dt>20, dt<100, limmag>=16.2))
plt.errorbar(dt[choose], f[choose], ef[choose], fmt='o')

dat = dat[choose]
choose = filt=='c'
errorbar(dat['###MJD'][choose], dat['uJy'][choose], dat['duJy'][choose], c='cyan', fmt='o')
choose = filt=='o'
errorbar(dat['###MJD'][choose], dat['uJy'][choose], dat['duJy'][choose], c='orange', fmt='o')

sig = dat['uJy']/dat['duJy']
dat['mag5sig Sky'].values

# 18cow
dat = pd.read_fwf("../data/opt/at2018cow_photometry_table.dat")
t = dat['MJD'].values
choose = t > 58285
t = t[choose].astype(float)
mag = dat['ABMag'][choose].values.astype(float)
emag = dat['Emag'][choose].values.astype(float)
dt = t-t[0]
dt_hr = dt*24
dt_min = dt_hr*60
# bin into ten-min blocks
dt_ten_min = np.unique((dt_min/10).astype(int))
# bin by hour
dt_hr_int = np.unique((dt_hr.astype(int)))
