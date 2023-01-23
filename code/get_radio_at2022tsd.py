""" Get the radio observations """

import pandas as pd

def get_radio():
    dd = "/Users/annaho/Dropbox/astro/papers/papers_active/AT2022tsd/data/radio"
    vla = pd.read_csv(dd+"/vla.txt")

    noema = pd.read_csv(dd+"/noema.txt")
    # Convert to mJy
    noema['Flux'] = noema['Flux']/1000
    noema['eFlux'] = noema['eFlux']/1000

    alma = pd.read_csv(dd+"/alma.txt")

    radio = pd.concat([vla,noema,alma],axis=0,ignore_index=True).sort_values('Date', ignore_index=True)
    return radio

