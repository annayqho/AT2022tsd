""" Calculate the duration of the optical LC """
import numpy as np
import matplotlib.pyplot as plt
from get_opt import *
import vals


def bin_lc(x,y,ey,bin_size):
    """ Bin a light curve into windows of size bin_size
    Return the binned light curve
    """
    x_binned = []
    y_binned = []
    ey_binned = []

    for ii,x_val in enumerate(x):
        choose = np.abs(x-x_val)<bin_size
        if sum(choose)==1:
            x_binned.append(x_val)
            y_binned.append(y[choose][0])
            ey_binned.append(ey[choose][0])
        elif sum(choose)>1:
            mean,wsum = np.average(
                y[choose], weights=1/ey[choose]**2, returned=True)
            efmean = np.sqrt(1/wsum)
            x_binned.append(np.average(x[choose]))
            y_binned.append(mean)
            ey_binned.append(efmean)

    x_binned = np.array(x_binned)
    y_binned = np.array(y_binned)
    ey_binned = np.array(ey_binned)

    return x_binned,y_binned,ey_binned


def single_band_flux_fit(ax,t,y,ey,zval,units='flux'):
    """ Do the fit for a single observing band """
    # Plot LC
    peak_ind = np.argmax(y)
    if units=='mag':
        peak_ind = np.argmin(y)
    jd_peak = t[peak_ind]
    f_peak = y[peak_ind]
    dt = t-jd_peak # in observer-frame
     
    # Do all calculations and plotting in the rest frame
    dt = dt/(1+zval)

    ax.errorbar(dt,y,ey,fmt='o',c='k')

    if units=='flux':
        halfmax = f_peak/2
    elif units=='mag':
        halfmax = f_peak+0.75
    print(halfmax)
    ax.axhline(y=halfmax, lw=0.5)

    # Initialize the values
    trise_str = ''
    tfade_str = ''

    # Check if the rise is resolved
    last_point = np.logical_and(y<halfmax, dt<=0)
    first_point = np.logical_and(y>halfmax, dt<=0)
    # if it's within 2-sigma of being zero...
    if np.logical_or(y[last_point][-1]-2*ey[last_point][-1]<0,ey[last_point][-1]==0):
        # if the closest point to half-max is a non-detection,
        # the rise is just the time between the detections
        lastdt = dt[last_point][-1]
        firstdt = dt[first_point][0]
        if firstdt==0:
            trise_str = '<%s' %(np.round(-lastdt,2))
        else:
            trise_str = '%s-%s' %(np.round(-firstdt,2),np.round(-lastdt,2))

    # Check if the fade is resolved
    last_point = np.logical_and(y<halfmax, dt>=0)
    if sum(last_point)==0:
        lastdt = dt[-1]
        tfade_str = '>%s' %(np.round(lastdt,2))

    # Grid for interpolating
    grid = np.linspace(min(dt), max(dt), 10000)

    # Monte Carlo the light curve
    nsim = 500
    trise = np.zeros(nsim)-99
    tfade = np.zeros(nsim)-99
    fpeak = np.zeros(nsim)-99
    tpeak = np.zeros(nsim)-99

    ysamples = np.zeros((nsim,len(t)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    # For each light curve draw, calculate rise time, fade time, and peak mag
    for ii in np.arange(nsim):
        ysample = ysamples[ii]

        # Identify the peak
        peak_ind = np.argmax(ysample)
        if units=='mag':
            peak_ind = np.argmin(ysample)

        jd_peak = t[peak_ind]
        tpeak[ii] = jd_peak
        f_peak = ysample[peak_ind]
        fpeak[ii] = f_peak

        lim = f_peak/2
        if units=='mag':
            lim = f_peak+0.75

        # Create an interpolated light curve
        yvals = np.interp(grid,dt,ysample)
        ax.plot(grid,yvals,lw=0.1,alpha=0.1,c='k')

        # Calculate the rise time
        reg = grid[np.logical_and(yvals<lim,grid<0)]
        if units=='mag':
            reg = grid[np.logical_and(yvals>lim,grid<0)]
        if len(reg)>0:
            trise[ii] = (np.abs(max(reg)))

            # Calculate the fade time
            reg = grid[np.logical_and(yvals<lim,grid>0)]
            if units=='mag':
                reg = grid[np.logical_and(yvals>lim,grid>0)]
            if len(reg)>0:
                tfade[ii] = (min(reg))
            if len(reg)==0:
                tfade[ii] = np.max(yvals)

    rise = np.mean(trise[trise>-99])
    erise = np.std(trise[trise>-99])
    if trise_str == '':
        trise_str = '%s\pm%s' %(format(rise,'.2f'), format(erise,'.2f'))

    fade = np.mean(tfade[tfade>-99])
    efade = np.std(tfade[tfade>-99])
    if tfade_str== '':
        tfade_str = '%s\pm%s' %(format(fade,'.2f'), format(efade,'.2f'))

    peak = np.mean(fpeak)
    epeak = np.std(fpeak)

    if units=='flux':
        mpeak = -2.5*np.log10(peak/1E6)*8.90
        empeak = np.abs(-1.0855*epeak/peak)

    ax.axvline(x=fade,lw=0.5)
    ax.axvline(x=-rise,lw=0.5)
    ax.axhline(y=f_peak,lw=0.5)
    ax.set_title("rise: %s; fade: %s" %(trise_str, tfade_str))
    ax.set_xlim(-30,50)
    ax.set_xlabel("Days")
    ax.set_ylabel("Flux")

    if units=='mag':
        plt.gca().invert_yaxis()

    peak_jd = np.mean(tpeak)
    return trise_str,tfade_str


if __name__=="__main__":
    jd,exp,filt,mag,mag_extcorr,emag,f,ef,f_extcorr,ef_extcorr = get_ipac()
    choose = filt=='g'
    fig,ax = plt.subplots(1,1,figsize=(3,3)) 
    single_band_flux_fit(
            ax,jd[choose],f[choose],ef[choose],vals.z,units='flux')
    plt.show()
