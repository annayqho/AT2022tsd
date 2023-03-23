import numpy as np
from astropy.cosmology import Planck18
from astropy.coordinates import SkyCoord
import extinction
from dustmaps.sfd import SFDQuery
from dustmaps.planck import PlanckQuery
from dustmaps.config import config
config['data_dir'] = '/Users/annaho/Dropbox/astro/tools/scanning/sfddata-master/'


# Important values
ra = 50.045308 
dec = 8.748872
t0 = 2459829.9731713 # first ZTF detection (from IPAC forced phot)
last_nondet = 2459826.9464 # last ZTF non-detection
t_flare_onset = 2459856.9 # onset of flaring
z = 0.2567 # redshift
dm = Planck18.distmod(z=z).value
dL_cm = Planck18.luminosity_distance(z=z).cgs.value
dL_mpc = Planck18.luminosity_distance(z=z).value
# I got these colors from https://www.visualisingdata.com/2019/08/five-ways-to-design-for-red-green-colour-blindness/
rc = '#DB4325'
gc = '#57C4AD'
ic = '#EDA247'
uc = '#E6E1BC'
wc = 'grey'

# Get extinction
def get_extinction(wave):
    # wavelength in AA
    coords = SkyCoord(ra,dec,unit='deg')
    sfd = SFDQuery()
    ebv = sfd(coords)
    a_v = ebv*2.742 # conversion for V-band (Schlafly & Finkbeiner 2011)
    ext = extinction.fitzpatrick99(np.array([wave]),a_v,r_v=3.1,unit='aa')
    return ext[0]

ext = {}
ext['u'] = 1.162 # Schlafly & Finkbeiner (2011)
ext['g'] = 0.906 # Schlafly & Finkbeiner (2011)
ext['r'] = 0.627 # Schlafly & Finkbeiner (2011)
ext['R'] = 0.627 # Schlafly & Finkbeiner (2011) # Should double check this
ext['i'] = 0.466 # Schlafly & Finkbeiner (2011)
ztf_pivot = {}
ztf_pivot['g'] = 4758.76
ztf_pivot['r'] = 6389.72
ztf_pivot['i'] = 7927.51
sdss_pivot = {}
sdss_pivot['g'] = 4702.50
sdss_pivot['r'] = 6175.58
sdss_pivot['i'] = 7489.98
keck_leff = {}
keck_leff['u'] = 3450.0
keck_leff['g'] = 4706.0
keck_leff['i'] = 7599.0
ps1_leff = {}
ps1_leff['i'] = 7520.0
ps1_leff['w'] = 6080.0
ps1_leff['z'] = 8660.0
