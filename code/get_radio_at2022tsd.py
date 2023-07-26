""" Get the radio observations """

import pandas as pd
import numpy as np
from astropy.time import Time
import vals

def get_radio():
    dd = "../../data/radio"
    sma = pd.read_csv(dd+"/sma.txt")
    vla = pd.read_csv(dd+"/vla.txt")
    noema = pd.read_csv(dd+"/noema.txt")

    # Add :00 to the end of each time
    dates = vla['Date'].values
    vla['Date'] = np.array([d+":00" for d in dates])

    # Add :00 to the end of each time
    dates = noema['Date'].values
    noema['Date'] = np.array([d+":00" for d in dates])

    # Add :00 to the end of each time
    dates = sma['Date'].values
    sma['Date'] = np.array([d+":00" for d in dates])

    # Convert to mJy
    noema['Flux'] = noema['Flux']/1000
    noema['eFlux'] = noema['eFlux']/1000

    alma = pd.read_csv(dd+"/alma.txt")

    gmrt = pd.read_csv(dd+"/gmrt.txt")

    radio = pd.concat([sma,gmrt,vla,noema,alma],axis=0,ignore_index=True).sort_values('Date', ignore_index=True)
    return radio


def get_data():
    """ Get the data for plotting """
    # Load in the data
    dat = get_radio()
    dt = (Time(dat['Date'].values.astype(str),format='isot').jd-vals.t0)/(1+vals.z)
    dat['dt'] = dt

    # Need to add systematic uncertainties to the RMS noise
    # Add 10% to ALMA data
    flux = dat['Flux'].values
    eflux = dat['eFlux'].values
    choose = dat['Tel'].values=='ALMA'
    eflux[choose] = np.sqrt((dat['eFlux'][choose].values)**2+(0.1*flux[choose])**2)

    # Add 5% to Ku band
    choose = np.logical_and(dat['Tel'].values=='VLA', dat['Freq_Obs'].values==15)
    eflux[choose] = np.sqrt((dat['eFlux'][choose].values)**2+(0.05*flux[choose])**2)

    # Add 15% to the three higher bands, for detections
    choose = np.logical_and.reduce((
            dat['Tel'].values=='VLA',dat['Freq_Obs'].values>15,
            dat['Flux']<99))
    eflux[choose] = np.sqrt(
            (dat['eFlux'][choose].values)**2+(0.15*flux[choose])**2)

    dat['eFlux'] = eflux
    return dat
