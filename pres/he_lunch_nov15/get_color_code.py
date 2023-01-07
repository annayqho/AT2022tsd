""" Get the color code for various plots """

import numpy as np
import pandas as pd


def get_cc():
    b = pd.read_csv("basic_info.csv")
    cl = b['Class'].values

    # Marker styles for different classes
    shape = {}
    shape['Unknown'] = 'o'
    shape['featureless'] = 'o'
    shape['AT2018cow-like'] = 'D'
    shape['IIb'] = 's'
    shape['IIn'] = '^'
    shape['Ic-BL'] = 'X'
    shape['Ic'] = 'D'
    shape['Ibn'] = 'P'
    shape['Ib'] = '*'
    shape['II'] = 'p'
    shape['IIn/Ibn?'] = 'o'
    shape['IIn/Ibn'] = 'o'
    shape['Icn'] = 'D'

    ec = {}
    ec['featureless'] = 'k'
    ec['Unknown'] = 'lightgrey'
    ec['AT2018cow-like'] = 'Crimson'
    ec['IIb'] = '#5cd061'
    ec['II'] = '#9ec480'
    ec['Ic'] = '#208eb7'
    ec['IIn'] = '#8a238b'
    ec['Ic-BL'] = '#9ecbf4'
    ec['Ibn'] = '#0a4f4e'
    ec['Ib'] = '#4ad9e1'
    ec['IIn/Ibn?'] = '#2c457d'
    ec['IIn/Ibn'] = '#2c457d'
    ec['Icn'] = 'k'

    msize = {}
    msize['Unknown'] = 5
    msize['featureless'] = 5
    msize['AT2018cow-like'] = 10
    msize['IIb'] = 10
    msize['II'] = 10
    msize['IIn'] = 10
    msize['Ic-BL'] = 10
    msize['Ic'] = 9
    msize['Ibn'] = 10
    msize['Ib'] = 15
    msize['IIn/Ibn?'] = 10
    msize['IIn/Ibn'] = 10
    msize['Icn'] = 10

    # the only one where it's different is featureless objects
    fc = {}
    for val in np.unique(cl):
        if val=='featureless':
            print("yes")
            fc[val] = 'lightgrey'
        else:
            fc[val] = ec[val]
    fc['Icn'] = 'white'

    return ec,fc,msize,shape


def get_cc_lit():
    """ Get the color code for the literature sample """
    cl = np.array(['Ic','Ibn','KSN','Ic-BL','IIn/Ibn?','PS1','DES','SNLS','IIn'])

    ec_ztf, fc_ztf, msize_ztf, shape_ztf = get_cc()

    # Marker styles for different classes
    shape = {}
    shape['Ibn'] = shape_ztf['Ibn']
    shape['KSN'] = '*'
    shape['Ic-BL'] = shape_ztf['Ic-BL']
    shape['IIn/Ibn?'] = shape_ztf['IIn/Ibn?']
    shape['PS1'] = 'o'
    shape['DES'] = 's'
    shape['SNLS'] = 'D'
    shape['IIn'] = shape_ztf['IIn']
    shape['Ic'] = 'D'

    ec = {}
    ec['IIn'] = ec_ztf['IIn']
    ec['Ic-BL'] = ec_ztf['Ic-BL']
    ec['Ibn'] = ec_ztf['Ibn']
    ec['IIn/Ibn?'] = ec_ztf['IIn/Ibn?']
    ec['KSN'] = 'k'
    ec['PS1'] = 'k'
    ec['DES'] = 'k'
    ec['SNLS'] = 'k'
    ec['Ic'] = '#208eb7'

    msize = {}
    msize['IIn'] = msize_ztf['IIn']
    msize['Ic-BL'] = msize_ztf['Ic-BL']
    msize['Ibn'] = msize_ztf['Ibn']
    msize['IIn/Ibn?'] = msize_ztf['IIn/Ibn?']
    msize['KSN'] = 10
    msize['DES'] = 5
    msize['PS1'] = 5
    msize['SNLS'] = 5
    msize['Ic'] = 9

    fc = {}
    for val in np.unique(cl):
        fc[val] = ec[val]

    return ec,fc,msize,shape 
