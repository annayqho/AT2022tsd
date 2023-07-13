""" read the Kann light curves """

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


def load_kann_lc():
    inputf = open("./SAMPLE!!! with Galactic Extinction.csv", "r")
    raw = np.array(inputf.readlines())
    inputf.close()

    stripped = np.array([val.rstrip() for val in raw])
    hascommas = np.array([',,' in val for val in stripped])
    while sum(hascommas) > 0:
        stripped = np.array([str.replace(val, ',,', ',') for val in stripped])
        hascommas = np.array([',,' in val for val in stripped])

    filt = np.array([list(filter(None, val.split(','))) for val in stripped], 
                    dtype='object')
    lengths = np.array([len(val) for val in filt])
    keep = filt[lengths > 0]

    # OK, first: we want everything that's a long GRB, not a short GRB.
    # Stopping criterion: the line with the string 'short'
    isshort = np.array(['Short' in val[0] for val in keep])
    stop = np.where(isshort)[0][0]

    long_grbs = keep[0:stop]
    # delete all rows that contain '---' or 'Swift'
    bad = np.array(['Swift' in val[0] or '--' in val[0] for val in long_grbs])
    long_grbs = long_grbs[~bad]

    lengths = np.array([len(val) for val in long_grbs])

    # Now, start storing the individual light curves. Will be a pain.
    # Make a dictionary

    names = []
    redshifts = []
    name_inds = []

    for ii,row in enumerate(long_grbs):
        if np.logical_and(len(row) == 2, len(row[0]) < 10):
            names.append(row[0])
            redshifts.append(row[1])
            name_inds.append(ii)

    lc = dict.fromkeys(names)

    names = np.array(names)
    name_inds = np.array(name_inds)

    starts = (name_inds + 1)
    ends = name_inds[1:]
    ends = np.append(ends, len(long_grbs))

    for ii,name in enumerate(names):
        try:
            data = long_grbs[starts[ii]:ends[ii]]
            tuples = [tuple(val) for val in data] # necessary for np array
            tuples = np.array(tuples)
            t = tuples[:,0]
            mag = tuples[:,1]
            emag = tuples[:,2]
            lc[name] = {'z':redshifts[ii],'t':t,'mag':mag,'emag':emag}
        except:
            lc[name] = {'z':0,'t':[0],'mag':[0],'emag':[0]}
    
    return lc


def plot_lc():
    lc = load_lc()
    yes = []
    for key,value in lc.items():
        x = value['t'].astype(float)*24
        y = value['mag'].astype(float)
        if np.interp(1, x, y) <= 19:
            col = 'red'
            yes.append(True)
        else:
            col = 'k'
            yes.append(False)
        plt.plot(x, y, c=col, alpha=0.1) 
    plt.gca().invert_yaxis()
    plt.xscale('log')
    # 3 hours in days, 6 hours in days
    plt.axvline(x=3.0/24, c='k', lw=1)
    plt.axvline(x=6.0/24, c='k', lw=1)
    plt.show()


def choose_bright():
    """
    choose all GRBs bright enough at t = 1 day to be seen with ZTF
    """
    lc = load_lc()
    bursts = {}
    lim_mag = 20.5
    for key,value in lc.items():
        f = interp1d(value['t'].astype(float),value['mag'].astype(float))
        try:
            lc_day = f(1)
            if lc_day < 20.5:
                bursts[key] = value
        except:
            pass
    return bursts
