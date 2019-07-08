import matplotlib.pyplot as plt
from astroquery.mast import Observations
from mastURL import gen_url, dload_url


ids = ['126602751', '12421477', '382523882']
for this_id in ids:
    obsTable = Observations.query_criteria(target_name=this_id, obs_collection="HLSP", filters="TESS",
                                                     t_exptime=[1799, 1801])
    sector = obsTable['sequence_number']
    print(sector)
