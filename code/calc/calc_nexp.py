""" Calculate the number of exposures obtained on how many different nights
following the original IMACS flare detection """

import pandas as pd
from get_opt import *

# Note: this doesn't include LAST
dat = get_full_opt()

mjd_imacs = max(dat['mjdstart'][dat['#instrument']=='Magellan/IMACS'])

choose = np.logical_and(dat['mjdstart']>mjd_imacs, dat['#instrument']!='ZTF')
dat = dat[choose]

# Total exposure time in seconds
exptot = sum(dat['exp'].values.astype(float)) 

# Total number of nights
nnights = len(np.unique(dat['mjdstart'].values.astype(int)))

# Now, add the LAST information
datlast = get_last()
exptot = exptot + 10.9*60*60 # 10.9 hours
nnights = len(np.unique(np.hstack((np.unique(dat['mjdstart'].values.astype(int)), np.unique(datlast[0].astype(int))))))
# Number of different telescopes
ntel = len(np.unique(dat['#instrument'].values)) # one of these is PS1, and add LAST...so they cancel out
