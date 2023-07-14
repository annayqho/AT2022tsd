""" Print details on each of the flare searches """

import sys
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/code")
import vals
from get_opt import *


def get_tab():
    # Retrieve the data
    tab = get_full_opt()

    # Starting place in the table: the start of the Magellan/IMACS sequence
    tstart = min(tab['mjdstart'][tab['#instrument']=='Magellan/IMACS'])
    tab = tab[tab['mjdstart']>tstart].reset_index()

    # Don't include ZTF
    tab = tab[tab['#instrument']!='ZTF'].reset_index()
    return tab


def print_telescope_data(tab):
    """ Instead of printing each sequence, print a summary for each telescope """
    # get unique telescopes, in order
    tel = list(dict.fromkeys(tab['#instrument'].values))

    # for each telescope, get total number of exposures and exposure time
    for t in tel:
        choose = tab['#instrument'].values==t
        exptot = int(sum(tab[choose]['exp'].values.astype(float))/60)
        filts = ''.join(np.unique(tab[choose]['flt']))
        nexp = sum(tab[choose]['nobs'])
        depth = np.round(np.median(tab[choose]['maglim']), 1)
        out = np.array([t, "$"+filts+"$", nexp, exptot, depth])
        print(' & '.join(out) + " \\\\")


def print_all_sequences(tab):
    """ Print all sequences """
    # I think we want to sort the table first by Instrument, then by date
    # tab = tab.sort_values(by=['#instrument', 'mjdstart'])
    # But then if you do that, you get GIT first, which is kind of weird

    # Starting conditions (already have tstart above)
    inst = 'Magellan/IMACS'
    filt = ['g']
    nexp = 0
    texp = 180
    maglims = []
    hasflare = False

    for index,row in tab.iterrows():
        # For each row, check whether it has the same instrument
        # and if the times roughly match up
        # If so, continue to add information
        is_next = row['mjdstart'] < tstart + texp/86400 + 15/60/24 # 15 min buffer
        if np.logical_and(row['#instrument']==inst, is_next):
            tstart = row['mjdstart']
            nexp += row['nobs']
            maglims.append(row['maglim'])
            filt.append(row['flt'])
            if row['isflare']:
                hasflare = True
        # If you're done with that sequence, print its information
        else:
            maglims = np.array(maglims)
            filt = np.array(filt)
            if len(np.unique(filt))==1:
                filtstr = "$" + filt[0] + "$"
            else:
                filtstr = "$" + ''.join(filt) + "$"
            toprint = [str(np.round(tstart,4)),inst,filtstr,str(nexp),str(texp), 
                       str(np.round(np.median(maglims),2)), str(hasflare)]
            print(" & ".join(toprint) + "\\\\")
            # Reset values
            tstart = row['mjdstart']
            inst = row['#instrument']
            filt = [row['flt']]
            nexp = row['nobs']
            texp = row['exp']
            maglims = []
            hasflare = False


if __name__=="__main__":
    tab = get_tab()
    print_telescope_data(tab)


    

