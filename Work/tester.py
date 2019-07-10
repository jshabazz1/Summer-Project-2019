from flatwrm import flatwrm
from matplotlib import pyplot as plt
from Time_and_flux import gettimeflux_1800
from astroquery.mast import Observations
from get_IDS import get_IDS
from gridPlots import flare_focus
import numpy as np




hardcopy = False
outfile = ""
debug = False
noplot = False
fit_events = True
magnitude = False
flarepoints = 3
sigma = 5
period=0.
degree=0
fwhm = 0.

#get_IDS(10)

ids = get_IDS(100)#['167602025']#, '167602025', '167814740']
for this_id in ids:
    try:
        sectors_search = Observations.query_criteria(target_name=this_id, obs_collection="HLSP", filters="TESS",
                                                     t_exptime=[1799, 1801])
        if  len(sectors_search):
            sector = sectors_search['sequence_number']

        time = gettimeflux_1800(this_id, str(sector[0]))[0]
        flux = gettimeflux_1800(this_id, str(sector[0]))[1]

        if period == 0:
            period = flatwrm.FindPeriod(time, flux, maxper=10,debug=debug)
        window_var = 1.5

        if degree == 0:
            degree = flatwrm.SelectDegree(time, flux, period * window_var, debug=debug)
        istart, istop = flatwrm.FindFlares(time, flux, period,
                               returnbinary=False,
                               N3=flarepoints,
                               degree=degree,
                               detection_sigma=sigma,
                               debug=debug)

        if hardcopy or outfile != "":
            if outfile == "":
                outfile = this_id + ".flare"
            if debug:
                print("Saving output into:", outfile)

        flatwrm.GenerateOutput(time, flux, istart, istop, period, \
                       fit_events=fit_events, degree=degree, debug=debug, outputfile=outfile)
        outfile = ""
        period = 0.
        degree = 0
        fwhm = 0.


        numFlares = len(istart)

    except:
        continue




    #plt.scatter(time, flux)
    #plt.show()

    #plt.scatter(time, flux, c='m')
    #plt.subplots(figsize = (12,7))
    from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

    figs, ax = plt.subplots()
    ax.scatter(time, flux)
    ax.scatter(time[istart], flux[istart], c = 'm')
    ax.scatter(time[istop], flux[istop], c = 'm')
    figs.suptitle("TIC " + this_id + ": Flare candidates")
    for t0, t1 in zip(istart, istop):
        ax.scatter(time[t0:t1+1], flux[t0:t1+1], c='m')
    plt.show()
#    plt.show()
#    axins = zoomed_inset_axes(ax, 2.5, loc=2)
#   axins.plot(time, flux, 'ro')
#    x1, x2, y1, y2 = time[istart][0], time[istop][0], flux[istart][0], flux[istop][0]
#    axins.set_xlim(x1-1.0, x2+1.0)
#    axins.set_ylim(y1-1, y2)
#    plt.show()
#
#
    fig, axs = plt.subplots(2, 2, figsize = (12,7))
    for istart, istop in zip(istart, istop):
        axs[0, 0].scatter(time, flux, c='k')
        axs[0, 0].scatter(time[istart:istop+1], flux[istart:istop+1])
        axs[0, 0].set_xlim(time[istart]-1, time[istop]+1)
        axs[0, 0].set_title('Flare candidate 1')

        axs[0, 1].scatter(time, flux, c='k')
        axs[0, 1].scatter(time[istart:istop+1], flux[istart:istop+1])
        axs[0, 1].set_xlim(time[istart]-1, time[istop]+1)
        axs[0, 1].set_title('Flare candidate 2')

        axs[1, 0].scatter(time[istart:istop+1], flux[istart:istop+1])
        #axs[1, 0].scatter(time, flux, c='k')
        #axs[1, 0].xlim(time[istart[2]]-2, time[istop[2]]+2)
        axs[1, 0].set_title('Flare candidate 3')

        axs[1, 1].scatter(time[istart:istop+1], flux[istart:istop+1])
        #axs[1, 1].scatter(time, flux, c='k')
        #axs[1, 1].xlim(time[istart[3]]-2, time[istop[3]]+2)
        axs[1, 1].set_title('Flare candidate 4')
    plt.show()













