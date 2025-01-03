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
    assert np.all(np.isclose(anc.rh, 71.3460312, atol=1e-7))
    assert np.all(np.isclose(anc.ws, 3.78513718, atol=1e-7))
    assert np.all(np.isclose(anc.wd, 36.5364418, atol=1e-7))
    assert np.all(np.isclose(anc.pr, 1018.85895, atol=1e-7))

    # precipitable water (water vapor)
    assert np.all(np.isclose(anc.wv, 2.63046503, atol=1e-7))

    # ozone
    # assert np.all(np.isclose(anc.oz, 0.268703669, atol=1e-7))

    # no2 and fraction
    # assert np.all(np.isclose(anc.no2_tropo, 3.50783159e15, atol=1e-7))
    # assert np.all(np.isclose(anc.no2_strat, 2.74156191e15, atol=1e-7))
    # assert np.all(np.isclose(anc.no2_frac, 0.53581202, atol=1e-7))
