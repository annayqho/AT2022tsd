""" Print the list of authors, affiliations, and acknowledgements """

import numpy as np
import pandas as pd


def print_authors():
    """ Print the list of authors and corresponding list of affiliations """
    # Load data
    dat = pd.read_excel("../data/coauthors.xlsx")

    # Add column for surname
    dat['Surname'] = [val.split(' ')[-1] for val in dat['Name'].values]

    ### Put all the affiliations in the correct order in the array
    affs = []

    # First, go through the ordered co-authors
    ordered = dat[~np.isnan(dat['Order'].values)].sort_values('Order').reset_index()
    for index,row in ordered.iterrows():
        ack = row['Acknowledgements']
        if pd.isnull(ack)==False:
            print(ack)
        for i in [1,2,3]:
            affil = row['Affil%s'%i]
            affs.append(affil)

    # Next, go through the rest of the co-authors
    alphabetical = dat[np.isnan(dat['Order'].values)].sort_values('Surname').reset_index()
    for index,row in alphabetical.iterrows():
        ack = row['Acknowledgements']
        if pd.isnull(ack)==False:
            print(ack)
        for i in [1,2,3]:
            affil = row['Affil%s'%i]
            affs.append(affil)

    # Remove nans
    affs = np.array(affs).astype(str)
    affs = affs[affs!='nan']
    # Make unique array
    _,idx = np.unique(affs, return_index=True)
    affs_unique = affs[np.sort(idx)]
        
    # Now, go through authors in order and give exponent of affil #
    authorstr = ""
    for index,row in ordered.iterrows():
        name = row['Name']
        authorstr += "%s$^{"%(name)
        for i in [1,2,3]:
            affil = row['Affil%s'%i]
            if pd.isnull(affil)==False:
                ind = np.where(affs_unique==affil)[0][0]
                authorstr += "%s,"%(ind+1)
        authorstr = authorstr[:-1]
        authorstr += "}$, "

    for index,row in alphabetical.iterrows():
        name = row['Name']
        authorstr += "%s$^{"%(name)
        for i in [1,2,3]:
            affil = row['Affil%s'%i]
            if pd.isnull(affil)==False:
                ind = np.where(affs_unique==affil)[0][0]
                authorstr += "%s,"%(ind+1)
        authorstr = authorstr[:-1]
        authorstr += "}$, "
    authorstr = authorstr[:-2]
    # Have to remove triple backslashes by hand

    # Produce affiliation string
    for aff in affs_unique:
        print("\\item %s" %aff)


def print_acknowledgements():
    # Load data
    dat = pd.read_excel("../data/acknowledgements.xlsx")
    for index,row in dat.iterrows():
        for i in [1,2]:
            ack = row['Ack%s'%i]
            if pd.isnull(ack)==False:
                print(ack)
                print("") # \n doesn't work in latex
