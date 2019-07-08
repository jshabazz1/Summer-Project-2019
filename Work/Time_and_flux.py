import mastURL
from astropy.io import fits


def gettimeflux_1800(tic_id, sector):
    url = mastURL.gen_url(tic_id, sector)
    fits.getdata(url, ext=1)
    with fits.open(url, mode="readonly") as hdulist:
        time = hdulist[1].data['TIME']
        flux = hdulist[1].data['FLUX_RAW']
    return time, flux

def gettimeflux_120(tic_id, sector):
    manifest = mastURL.dload_url(tic_id, sector)
    url = manifest[1][0]
    fits.getdata(url, ext=1)
    with fits.open(url, mode="readonly") as hdulist:
        time = hdulist[1].data['TIME']
        flux = hdulist[1].data['PDCSAP_FLUX']
    return time, flux


