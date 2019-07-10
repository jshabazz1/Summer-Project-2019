from astroquery.mast import Catalogs


def get_IDS (number):

    catalogTIC = Catalogs.query_criteria(catalog = "Tic", Tmag = [12.0,12.5], Teff = [3500.0, 3550.0] , logg = [4.2,5.0])
    catalogTIC = catalogTIC[catalogTIC['dec'] < 0]
    return catalogTIC['ID'][:number]






