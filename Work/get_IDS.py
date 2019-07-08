from astroquery.mast import Catalogs


def get_IDS (number):

    catalogTIC = Catalogs.query_criteria(catalog = "Tic", Tmag = [11.5,12.0], Teff = [3450.0, 3500.0] , logg = [4.2,5.0])
    catalogTIC = catalogTIC[catalogTIC['dec'] < 0]
    return catalogTIC['ID'][:number]






