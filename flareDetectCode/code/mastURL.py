from astroquery.mast import Observations



def dload_url(tic_id, sector): #Dowloads URL through query for general 2-minute cadence data given TIC ID and Sector
     obsTable = Observations.query_criteria(
            obs_collection = "TESS",
            dataproduct_type = ["TIMESERIES"],
            target_name = tic_id,
            sequence_number=sector)
     products = Observations.get_product_list(obsTable)
     manifest = Observations.download_products(products, extension = ['fits'])
     return manifest


def gen_url(tic_id, sector, version = "v04"):  # Generates URL for TASOC data given TIC ID, Sector, Cadence, and Version
    tic_id_1 = tic_id.zfill(16)
    tic_id_2 = tic_id.zfill(11)
    cadence = '1800'
    sub_dir = "ffi"


    url = "http://archive.stsci.edu/hlsps/tasoc/s000" + str(sector) + "/" + sub_dir + "/" \
          + tic_id_1[0:4] + "/" + tic_id_1[4:8] + "/" + tic_id_1[8:12] + "/" + tic_id_1[12:16] + "/" \
          + "hlsp_tasoc_tess_" + sub_dir + "_tic" + tic_id_2 + "-s0" + str(sector) + "-c" + cadence + "_tess_" \
          + version + "_lc.fits"
    return url




