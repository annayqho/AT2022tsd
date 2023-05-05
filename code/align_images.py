""" Align images using SWARP """

import os
import vals

swarp_path = '/opt/local/bin/swarp'

dd = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/opt/LRIS/imaging/"
u1 = dd+"aeazfpb221229_0081r.fits" # has a flare
u2 = dd+"aeazfpb221229_0082r.fits"
u3 = dd+"aeazfpb221229_0083r.fits"
u4 = dd+"aeazfpb221229_0084r.fits"
u5 = dd+"aeazfpb221229_0085r.fits"
g = dd+"lrisb20230117_ZTF22abftjko_G.fits"
i = dd+"lrisr20230117_ZTF22abftjko_i.fits"
u = dd+"coadd_u.fits"

# centering
ra = vals.ra
dec = vals.dec
size = 2000

def coadd_uband():
    swarp_command = swarp_path \
            + " %s %s %s %s -c config.swarp -CENTER '" %(u2,u3,u4,u5)\
            + str(ra) + " " + str(dec) \
            + "' -SUBTRACT_BACK Y -RESAMPLE Y -COMBINE Y -IMAGE_SIZE '" \
            + str(size) + "," + str(size) + "'"
    os.system(swarp_command)


def align_rgb():
    swarp_command = swarp_path \
            + " %s %s %s -c config.swarp -CENTER '" %(g,i,u)\
            + str(ra) + " " + str(dec) \
            + "' -SUBTRACT_BACK Y -RESAMPLE Y -COMBINE N -IMAGE_SIZE '" \
            + str(size) + "," + str(size) + "'"
    os.system(swarp_command)


if __name__=="__main__":
    align_rgb()
