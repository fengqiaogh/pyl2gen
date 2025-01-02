import math

from oel_hdf4.libnav.nav import jd


def esdist(iyr, iday, msec):
    """
    This function computes the earth-sun distance in AU. It uses the model
    referenced in The Astronomical Almanac for 1984, Section S (Supplement).

    Arguments:

    iyr : int
         Year, four digits (i.e, 1993)
    iday : int
         Day of year (1-366)
    msec : int
         milliseconds of day

    Returns:
    rs : float
         Magnitude of the Sun vector (AU)
    """

    radeg = 57.29577951
    imon = 1

    # Compute floating point days since Jan 1.5, 2000
    # Note that the Julian day starts at noon on the specified date
    t = jd(iyr, imon, iday) - 2451545.0 + (msec / 1000.0 - 43200.0) / 86400.0

    # Compute mean anomaly
    gs = 357.52772 + 0.9856002831 * t

    # Compute solar distance (AU)
    res = (
        1.00014 - 0.01671 * math.cos(gs / radeg) - 0.00014 * math.cos(2.0 * gs / radeg)
    )

    return res
