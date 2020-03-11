from astroquery.mast import Catalogs
from astropy.table import Table


def get_IDS (number1, number2):

    catalogTIC = Catalogs.query_criteria(catalog = "Tic", Teff = [3450.0, 3500.0], logg=[4.2,5.0], Tmag=[11.5,12.0])
    catalogTIC = catalogTIC[catalogTIC['dec'] < 0]
    return catalogTIC['ID'][number1:number2]




# catalogTIC = Catalogs.query_criteria(catalog = "Tic", Teff = [3450.0, 3500.0], logg=[4.2,5.0], Tmag=[11.5,12.0])
# catalogTIC = catalogTIC[catalogTIC['dec'] < 0]
# print(catalogTIC)