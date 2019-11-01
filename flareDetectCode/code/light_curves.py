import matplotlib.pyplot as plt
from astroquery.mast import Observations
from scipy import signal
from Time_and_flux import gettimeflux_1800, gettimeflux_120




ids = ['415438368', '450061615']
cadence = '0120'
sector = '5'
for this_id in ids:
    if (cadence == '1800' ):
        obsTable = Observations.query_criteria(target_name=this_id, obs_collection="HLSP", filters="TESS",
                                               t_exptime=[1799, 1801])
        if sector == '':
            sector = obsTable['sequence_number']
        else:
            continue
        tess_bjds = gettimeflux_1800(this_id, str(sector[0]))[0]
        sap_fluxes = gettimeflux_1800(this_id, str(sector[0]))[1]
        med_flux_1 = signal.medfilt(sap_fluxes, kernel_size=65)
    else:
        tess_bjds = gettimeflux_120(this_id, sector)[0]
        sap_fluxes = gettimeflux_120(this_id, sector)[1]
        med_flux_1 = signal.medfilt(sap_fluxes, kernel_size=101)


    fig, ax = plt.subplots()

    ax.plot(tess_bjds, sap_fluxes, 'bo')#RAW blue


    fig.suptitle("TIC " + this_id + " Flare")
    ax.set_ylabel("Flux (e-/s)")
    ax.set_xlabel("Time (TBJD)")

    plt.subplots_adjust(left=0.15)
    plt.show()

    detrend_lc = (sap_fluxes/med_flux_1) - 1


    fig, ax = plt.subplots()

    ax.plot(tess_bjds, detrend_lc, 'ko')

    fig.suptitle("TIC " + this_id + " Flare")
    ax.set_ylabel("PDCSAP Flux (e-/s) *filtered*")
    ax.set_xlabel("Time (TBJD)")

    plt.subplots_adjust(left=0.15)
    plt.show()
