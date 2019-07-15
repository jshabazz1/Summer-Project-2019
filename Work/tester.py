from flatwrm import flatwrm
from matplotlib import pyplot as plt
from Time_and_flux import gettimeflux_1800
from astroquery.mast import Observations
from get_IDS import get_IDS
from astropy.time import Time





hardcopy = False
outfile = ""
debug = False
noplot = False
fit_events = True
magnitude = False
flarepoints = 3
sigma = 4
period=0.
degree=0
fwhm = 0.


ids = [





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



    if numFlares != 0 :

        """Displays a light curve with highlighted flare candidates"""
        figs, ax = plt.subplots()
        ax.scatter(time, flux)
        ax.scatter(time[istart], flux[istart], c = 'm')
        ax.scatter(time[istop], flux[istop], c = 'm')
        figs.suptitle("TIC " + this_id + ": Flare candidates")
        for t0, t1 in zip(istart, istop):
            ax.scatter(time[t0:t1+1], flux[t0:t1+1], c='m')
        plt.show()
        plt.savefig('/Users/jshabazz/Work/lightcurves/' + this_id + '.png', bbox_inches='tight')


        """Displays a grid of plots focused in on each individual flare in the light curve"""
        fig, axs = plt.subplots(2, 2, figsize = (11,7))
        fig.suptitle(f'TIC {this_id}')
        for i in range(numFlares):
            if i > 3:
                print(f'Number of flares detected: {numFlares}')
                break
            t2 = istart[i]
            t3 = istop[i]
            axs[i//2, i%2].scatter(time, flux, c='k')
            axs[i//2, i%2].scatter(time[t2:t3+1], flux[t2:t3+1])
            axs[i//2, i%2].set_xlim(time[t2]-0.5, time[t3]+0.5)
            axs[i//2, i%2].set_title(f'Flare candidate: {i+1}')
        plt.show()
        plt.savefig('/Users/jshabazz/Work/lightcurves/' + this_id + '_flares.png', bbox_inches='tight')
    else:
        continue















