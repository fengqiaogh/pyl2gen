from oel_hdf4.libnav.esdist import esdist
import numpy as np


def test_esdist():
    year = 2013
    day_of_year = 285
    milliseconds = 12875000
    res = esdist(year, day_of_year, milliseconds)
    assert np.isclose(res, 0.9980337, atol=1e-7)
