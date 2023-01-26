from astropy.cosmology import Planck18

# Important values
ra = 50.045308 
dec = 8.748872
t0 = 2459828.7257 # first GOTO detection
last_nondet = 2459826.9464 # last ZTF non-detection
z = 0.25666 # redshift
dm = Planck18.distmod(z=z).value
dL_cm = Planck18.luminosity_distance(z=z).cgs.value
dL_mpc = Planck18.luminosity_distance(z=z).value
rc = '#DB4325'
gc = '#57C4AD'
ext = {}
ext['g'] = 0.906 # Schlafly & Finkbeiner (2011)
ext['r'] = 0.627 # Schlafly & Finkbeiner (2011)
ext['i'] = 0.466 # Schlafly & Finkbeiner (2011)
ztf_pivot = {}
ztf_pivot['g'] = 4758.76
ztf_pivot['r'] = 6389.72
ztf_pivot['i'] = 7927.51
sdss_pivot = {}
sdss_pivot['g'] = 4702.50
sdss_pivot['r'] = 6175.58
sdss_pivot['i'] = 7489.98
