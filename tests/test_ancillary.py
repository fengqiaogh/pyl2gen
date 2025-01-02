from datetime import datetime
from pathlib import Path

import numpy as np

from ancillary import Ancillary


def test_set():
    lat = np.array([[31.0324459, 31.0324459], [31.0324459, 31.0324459]])
    lon = np.array([[122.401329, 122.401329], [122.401329, 122.401329]])
    dt = datetime(2013, 10, 12, 3, 16, 40)
    settings = {
        "met1": "share/common/met_climatology_v2014.h5",
        "ozone1": "share/common/ozone_climatology_v2014.h5",
        "no2file": "share/common/no2_climatology_v2013.h5",
        "gas_opt": -17,
        "OCDATAROOT": Path("share"),
    }
    anc = Ancillary()
    anc.set(lat, lon, dt, settings)
    anc.rh = None
    anc.ws = None
    anc.wd = None
    anc.pr = None
    anc.wv = None  # precipitable water (water vapor)
    anc.oz = None  # ozone

    # no2 and fraction
    anc.no2_tropo = None
    anc.no2_strat = None
    anc.no2_frac = None
    # assert np.isclose(results, 71.3460312, atol=1e-7)
