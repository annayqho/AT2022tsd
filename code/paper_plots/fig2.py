""" Combined Figure 2 

Required per Nature's editorial process """

import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "sans-serif"
from fig2_opt_lc import *

# Initialize
figwidth_mm = 183 # Nature standard
figwidth_in = (figwidth_mm/10)/2.54 # in inches
fig,axarr = plt.subplots(
        2,2,figsize=(figwidth_in,figwidth_in*(3/7)),
        sharey=True,gridspec_kw={'width_ratios': [1.5,3,2,2]})

plt.show()
