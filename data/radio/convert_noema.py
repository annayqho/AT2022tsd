import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code/")
from astropy.time import Time
import vals

# Start time
jd = []
rest_nu = []
obs_nu = []

with open("noema_raw.txt", "r") as inputf:
    lines = inputf.readlines()
    for line in lines:
        if 'UT' in line:
            dd = line[0:2]
            if line[3:6]=='OCT':
                mm = '10'
            elif line[3:6]=='NOV':
                mm = '11'
            elif line[3:6]=='DEC':
                mm = '12'
            else:
                print("bad month")
            tt = line[12:17]
            dat = '2022-%s-%sT%s' %(mm,dd,tt)
        elif 'GHz' in line:
            obs_freq = float(line.split(" ")[-2:][0])
        elif 'FLUX' in line:
            flux = line.split('=')[1].split('(')[0].replace(" ","")
            eflux = line.split('(')[1].split(')')[0].replace(" ","")
            print(dat+","+'NOEMA'+","+str(obs_freq)+","+str(flux)+","+str(eflux))
