""" Print the list of authors, affiliations, and acknowledgements """

import numpy as np
import pandas as pd

# Load data
dat = pd.read_excel("../data/coauthors.xlsx")

### Put all the affiliations in the correct order in the array
affs = []

# First, go through the ordered co-authors
ordered = dat[~np.isnan(dat['Order'].values)].sort_values('Order').reset_index()
for index,row in ordered.iterrows():
    for i in [1,2,3]:
        affil = dat['Affil%s'%i].values
        affs.append(affil)

# Next, go through the rest of the co-authors
alphabetical = dat[np.isnan(dat['Order'].values)].sort_values('Name').reset_index()
for index,row in ordered.iterrows():
    for i in [1,2,3]:
        affil = dat['Affil%s'%i].values
        affs.append(affil)

# Remove nans
aff = np.array(aff).astype(str)
aff = aff[aff!='nan']
# Make unique array
aff = np.unique(aff)
    
