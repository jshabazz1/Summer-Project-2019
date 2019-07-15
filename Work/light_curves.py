import matplotlib.pyplot as plt
from astroquery.mast import Observations
from scipy import signal
from Time_and_flux import gettimeflux_1800




ids = ['167602025', '167695269', '167695269', '167814740', '167814740' ]
cadence = '1800'
flare_times = [2458335.233, 2458336.636, 2458343.026, 2458344.888, 2458344.875]
for this_id in ids, bjd in flare_times:
    #for bjd in flare_times:
        tbjd = bjd - 2457000
        obsTable = Observations.query_criteria(target_name=this_id, obs_collection="HLSP", filters="TESS",
                                                         t_exptime=[1799, 1801])
        sector = obsTable['sequence_number']

        tess_bjds = gettimeflux_1800(this_id, str(sector[0]))[0]
        sap_fluxes = gettimeflux_1800(this_id, str(sector[0]))[1]

        if (cadence == '1800' ):
            med_flux_1 = signal.medfilt(sap_fluxes, kernel_size=65)
        else:
            med_flux_1 = signal.medfilt(sap_fluxes, kernel_size=101)


        t0 = tbjd #transit time

        fig, ax = plt.subplots()

        #ax.plot(tess_bjds, pdcsap_fluxes, 'ko')#BKG black
        #ax.plot(tess_bjds, med_flux_1, 'b-')
        ax.plot(tess_bjds, sap_fluxes, 'bo')#RAW blue



        ax.set_xlim(t0-1.0, t0+1.0)
        #ax.set_xlim(1335,1335.5)


        ax.axvline(x=t0, color = 'red')

        fig.suptitle("TIC " + this_id + " Flare")
        ax.set_ylabel("PDCSAP Flux (e-/s)")
        ax.set_xlabel("Time (TBJD)")

        plt.subplots_adjust(left=0.15)
        plt.show()

        detrend_lc = (sap_fluxes/med_flux_1) - 1


        fig, ax = plt.subplots()

        ax.plot(tess_bjds, detrend_lc, 'ko')


        ax.set_xlim(t0-0.5,t0+0.5)

        fig.suptitle("TIC " + this_id + " Flare")
        ax.set_ylabel("PDCSAP Flux (e-/s) *filtered*")
        ax.set_xlabel("Time (TBJD)")

        plt.subplots_adjust(left=0.15)
        plt.show()
        pass