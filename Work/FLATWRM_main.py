from flatwrm import flatwrm
from matplotlib import pyplot as plt
from Time_and_flux import gettimeflux_1800
import numpy as np



hardcopy = False
outfile = ""
debug = False
noplot = False
fit_events = True
magnitude = False
flarepoints = 2
sigma = 3
period=0.
degree=0
fwhm = 0.

ids = ['167695269', '167602025', '167814740']
for this_id in ids:
    time = gettimeflux_1800(this_id, '1')[0]
    flux = gettimeflux_1800(this_id, '1')[1]
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






    plt.scatter(time, flux)
    plt.show()

    plt.scatter(time, flux)
    for t0,t1 in zip(istart,istop):
        plt.scatter(time[t0:t1+1], flux[t0:t1+1], c='r')
    plt.show()

    for x in range(3):
        plt.scatter(time, flux)
        plt.scatter(time[istart][x], flux[istart][x], c='r')
        plt.scatter(time[istop][x], flux[istop][x], c='r')
        plt.xlim(time[istart][x]-2.0, time[istop][x]+2.0)
        plt.axvline((time[istart][x]+time[istop][x])/2, c='y')
        plt.show()


