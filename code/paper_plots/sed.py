""" Plot the X-ray to radio SED
Currently the only good date for this is 27-28d.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
import vals

# Initialize
fig,ax = plt.subplots(1,1,figsize=(4,3))

# Formatting
ax.set_xlabel(
        r"$\nu_\mathrm{obs}$",fontsize=11,
        fontname='sans-serif')

ax.set_ylabel(r"$L_\nu$", fontsize=11,
        fontname='sans-serif')

ax.tick_params(axis='both', fontsize=11)
plt.tight_layout()
plt.show()


