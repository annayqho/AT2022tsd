""" Check to see if any dwarf nova light curves resemble AT2022tsd """

from penquins import Kowalski
import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from helpers import *


if __name__=="__main__":
    # Get the names of the 182 CVs from BTS with mpeak > 18 mag
    cvs = pd.read_csv("../data/btscvs.txt")
    names = list(cvs['ZTFID'].values)

    # Log onto Kowalski
    s = logon()
    
    # Get IPAC command for each object
    for name in names:
        #print(name)
        ra,dec = get_pos(s, name)
        print(name,ra,dec)
        # JD start and end correspond to March 2018 and 1 Jan 2023
        #cmd = 'wget --http-user=ztffps --http-passwd=dontgocrazy! -O log.txt "https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?ra=%s&dec=%s&jdstart=2458178&jdend=2459945.5&email=annayqho@gmail.com&userpass=hasa654"' %(ra,dec)
        #print(cmd)
        # Print this to a file, then run using sh
