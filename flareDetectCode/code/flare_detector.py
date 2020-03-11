from astropy.table import Table
from flatwrm import flatwrm
from matplotlib import pyplot as plt
from Time_and_flux import gettimeflux_1800, get_catalog_data, get_header_data
from astroquery.mast import Observations
from get_IDS import get_IDS
import csv





hardcopy = False
outfile = "blurrrrrrrr.csv"
debug = False
noplot = False
fit_events = True
magnitude = False
flarepoints = 4
sigma = 8
period=0.
degree=0
fwhm = 0.


# ids = get_IDS(301, 1500)

guenther_flares = Table.read("/Users/jshabazz/flareDetectCode/textfiles/guenther_flares.csv",format='ascii')
ids = guenther_flares['tic_id'][0:760].astype(str)

#ids = ['141914082']
# with open("manifest_firstdraft.csv","w", newline='') as f:
#     fieldnames = ['TIC', 'Sector','Plot 1', 'Plot 2', 'Teff', 'Tmag', 'log(g)','CCD', 'Camera']
#     writer = csv.DictWriter(f, fieldnames=fieldnames)
#     writer.writeheader()
with open("test.csv", "w") as ofile, open("/Users/jshabazz/flareDetectCode/textfiles/manifest_guenther_1.csv","w", newline='') as f:
    # fieldnames = ['TIC_ID', 'Start_time', 'Stop_time', 'Start_flux', 'Stop_flux']
    # ofile.write(",".join(fieldnames)+'\n')
    fieldnames = ['TIC', 'Sector','Plot 1', 'Plot 2', 'Teff', 'Tmag', 'log(g)','CCD', 'Camera']
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for this_id in ids:
        try:
            #sectors_search = Observations.query_criteria(target_name=this_id, provenance_project='TASOC')
            sectors_search = Observations.query_criteria(target_name=this_id, obs_collection="HLSP", filters="TESS",
                                                         t_exptime=[1799, 1801])

            if  len(sectors_search):
                sector = sectors_search['sequence_number']


            time = gettimeflux_1800(this_id, str(sector[0]))[0]
            flux = gettimeflux_1800(this_id, str(sector[0]))[1]
            Tmag = get_catalog_data(this_id)[0]
            Teff = get_catalog_data(this_id)[1]
            logg = get_catalog_data(this_id)[2]
            CCD = get_header_data(this_id, str(sector[0]))[0]
            Camera = get_header_data(this_id, str(sector[0]))[1]

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
           
            flatwrm.GenerateOutput(time, flux, istart, istop, period, this_id, \
                           fit_events=fit_events, degree=degree, debug=debug, outputfile="")

            #save out put in readme text file


            period = 0.
            degree = 0
            fwhm = 0.


            numFlares = len(istart)

        except:
            continue

        #flatwrm.FitFlare(time, flux, istart, istop, period)

        if numFlares != 0 :
            

            """Displays a light curve with highlighted flare candidates"""
            figs, ax = plt.subplots()
            ax.scatter(time, flux)
            ax.scatter(time[istart], flux[istart], c = 'm')
            ax.scatter(time[istop], flux[istop], c = 'm')
            figs.suptitle("TIC " + this_id + ": Flare candidates")
            ax.set_xlabel('Time (TBJD)')
            ax.set_ylabel('Flux ')
            for t0, t1 in zip(istart, istop):
                ax.scatter(time[t0:t1+1], flux[t0:t1+1], c='m')
            plt.savefig('/Users/jshabazz/Work/guenther_flares_2/' + this_id + '.jpg', bbox_inches='tight')
            #plt.show()


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
                #axs[i//2, i%2].set_ylim(flux[t2], flux[t3]+50.0)
                axs[i//2, i%2].set_title(f'Flare candidate: {i+1}')
            plt.savefig('/Users/jshabazz/Work/guenther_flares_2/' + this_id + '_flares.jpg', bbox_inches='tight')
            #plt.show()
            writer.writerow({'TIC': this_id, 
                            'Plot 1': this_id + '.jpg', 
                            'Plot 2': this_id + '_flares.jpg', 
                            'Sector': str(sector[0]),
                            'Teff': str(Teff[0]),
                            'Tmag': str(Tmag[0]),
                            'log(g)': str(logg[0]),
                            'CCD': str(CCD),
                            'Camera': str(Camera)})
        else:
            continue















