from ztfquery import query
zquery = query.ZTFQuery()
jds = [2459829.9731713,2459829.9879167,2459840.9833912,2459856.9122454,2459857.8403241,2459881.6982407]


# Do the first images (with bright detections) as a test
zquery.load_metadata(radec=[50.045308, +8.748872],size=0.0001,
                     sql_query=f"obsjd>2459829 AND obsjd<2459830")

zquery.download_data("sciimg.fits")
zquery.download_data("scimrefdiffimg.fits.fz")

# Do the claimed i-band detection
# for some reason only returns one image
zquery.load_metadata(radec=[50.045308, +8.748872],size=0.0004,
                     sql_query=f"obsjd>2459854 AND obsjd<2459858")

# Do the claimed r-band detection at later times
zquery.load_metadata(radec=[50.045308, +8.748872],size=0.0004,
                     sql_query=f"obsjd>2459881.6 AND obsjd<2459881.7")

zquery.metatable['obsjd'].values
