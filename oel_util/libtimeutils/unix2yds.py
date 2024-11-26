import datetime
import math

from libgenutils.genutils import BAD_FLT, BAD_INT


def unix2yds(usec):
    """
    Converts secs since 1/1/70 to yr, day, secs of day.

    Args:
    usec : float
        seconds since 1/1/70

    Returns:
    tuple of (year, day, secs)
    year : int
        year
    day : int
        day of the year
    secs : float
        seconds of the day
    """
    if usec == BAD_FLT:
        year = BAD_INT
        day = BAD_INT
        secs = BAD_FLT
    else:
        utime = int(usec)

        # Use timezone-aware object for UTC datetime
        trec = datetime.datetime.fromtimestamp(utime, datetime.timezone.utc)

        year = trec.year
        day = trec.timetuple().tm_yday
        secs = trec.hour * 3600 + trec.minute * 60 + trec.second + math.fmod(usec, 1)
    return year, day, secs
