import eleanor
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

flares = Table.read("/Users/jshabazz/flareDetectCode/textfiles/guenther_flares.csv",format='ascii')
ids = flares['tic_id']
sector = flares['sector']
duration = flares['duration']

#%%

for i,this_id in enumerate(ids):
    if duration[i]>0.32:
        print(i)
        print(this_id)
        print(sector[i])
        star = eleanor.Source(tic = int(this_id) , sector = int(sector[i]))
        print('Found TIC {0} (Gaia {1}), with TESS magnitude {2}, RA {3}, and Dec {4}'.format(star.tic, star.gaia, star.tess_mag, star.coords[0], star.coords[1]))
        data = eleanor.TargetData(star, height=15, width=15, bkg_size=31, do_psf=True, do_pca=True)
        plt.figure(figsize = (15,5))

        q = data.quality == 0

        #plt.plot(data.time[q], data.raw_flux[q]/np.nanmedian(data.raw_flux[q])+0.06, 'k')
        #plt.plot(data.time[q], data.corr_flux[q]/np.nanmedian(data.corr_flux[q]) + 0.03, 'r')
        
        #plt.plot(data.time[q], data.psf_flux[q]/np.nanmedian(data.psf_flux[q]) - 0.02, 'b')
        plt.ylabel('Normalized Flux')
        plt.xlabel('Time [BJD - 2457000]')
        plt.title('TIC ' + this_id)
        
        flux_norm = np.average(data.pca_flux[q])
        std = np.std(data.pca_flux[q])

        outliers = []
        times = []
     
        

        for i in range(0, len(data.pca_flux[q])):
            if(data.pca_flux[q][i]>flux_norm + (3*std)):
                print('if statement')
                outliers.append(data.pca_flux[q][i])
                times.append(data.time[q][i])
        plt.plot(data.time[q], data.pca_flux[q]/np.nanmedian(data.pca_flux[q]), '-go')
        plt.plot(times, outliers/np.nanmedian(data.pca_flux[q]), '-ro')
        plt.show()
