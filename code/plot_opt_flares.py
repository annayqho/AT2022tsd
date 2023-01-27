""" Plot zoom-in of the optical flares to better understand what's going on

Might end up being paper figures later """

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import vals
from get_opt import *

tel,mjd,filt,mag,emag,limmag,flare = get_flares()

# Detections
is_det = ~np.logical_or(np.isnan(emag), mag>limmag)
dt = (Time(mjd, format='mjd').jd-vals.t0)/(1+vals.z)
plt.errorbar(dt[is_det], mag[is_det], emag[is_det], fmt='o', c='k')

# Upper limits
plt.scatter(
        dt[~is_det], mag[~is_det], marker='v', 
        edgecolor='k', facecolor='white')

plt.tight_layout()
plt.gca().invert_yaxis()
plt.show()



