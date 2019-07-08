from flatwrm import flatwrm
from matplotlib import pyplot as plt
from Time_and_flux import gettimeflux_1800
from astroquery.mast import Observations
from get_IDS import get_IDS




hardcopy = False
outfile = ""
debug = False
noplot = False
fit_events = True
magnitude = False
flarepoints = 4
sigma = 3
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



        #plt.scatter(time, flux)
        #plt.show()

        #plt.scatter(time, flux, c='m')
        plt.scatter(time[istart], flux[istart], c='b')
        plt.scatter(time[istop], flux[istop], c='b')
        plt.scatter(time, flux, c='k')
        plt.suptitle("TIC " + this_id + " Sector "+sector +": Flare candidates")
        plt.show
        for t0,t1 in zip(istart,istop):
            plt.scatter(time[t0:t1+1], flux[t0:t1+1], c='r')
        plt.show()

        #fig, axs = plt.subplots(2, numFlares)
        #for x in range(numFlares):
            #plt.scatter(time, flux)
            #for i in range(2):
                #axs[i].scatter(time[istart][x], flux[istart][x], c='r')
                #axs[i].scatter(time[istop][x], flux[istop][x], c='r')
                #plt[i].xlim(time[istart][x]-2.0, time[istop][x]+2.0)
                #plt[i].axvline((time[istart][x]+time[istop][x])/2, c='y')
                #plt[i].show()
    except:
        pass


