import numpy as np

from .gha2000 import gha2000
from .l_sun import sun2000


def sunangs(iyr, iday, gmt, xlon, ylat):

    pi = np.pi
    radeg = 180 / pi
    suni = np.zeros(3)
    sung = np.zeros(3)
    npix = xlon.size
    suna = np.zeros(npix)
    # up = np.zeros(3)
    up = np.zeros(shape=(npix, 3))

    # no = np.zeros(3)
    no = np.zeros(shape=(npix, 3))

    # ea = np.array([0.0, 0.0, 0.0])
    ea = np.zeros(shape=(npix, 3))

    sec = gmt * 3600.0
    suni, rs = sun2000(iyr, iday, sec)

    day = iday + sec / 86400.0
    gha = gha2000(iyr, day)
    ghar = gha / radeg

    sung[0] = suni[0] * np.cos(ghar) + suni[1] * np.sin(ghar)
    sung[1] = suni[1] * np.cos(ghar) - suni[0] * np.sin(ghar)
    sung[2] = suni[2]

    rlon = xlon / radeg
    rlat = ylat / radeg
    cosy = np.cos(rlat)
    siny = np.sin(rlat)
    cosx = np.cos(rlon)
    sinx = np.sin(rlon)

    up[:, 0] = cosy * cosx
    up[:, 1] = cosy * sinx
    up[:, 2] = siny

    upxy = np.sqrt(up[:, 0] * up[:, 0] + up[:, 1] * up[:, 1])
    ea[:, 0] = -up[:, 1] / upxy
    ea[:, 1] = up[:, 0] / upxy
    no[:, 0] = up[:, 1] * ea[:, 2] - up[:, 2] * ea[:, 1]
    no[:, 1] = up[:, 2] * ea[:, 0] - up[:, 0] * ea[:, 2]
    no[:, 2] = up[:, 0] * ea[:, 1] - up[:, 1] * ea[:, 0]

    sunv = np.dot(sung, up.T)

    sunn = np.dot(sung, no.T)
    sune = np.dot(sung, ea.T)
    sunz = radeg * np.arctan2(np.sqrt(sunn * sunn + sune * sune), sunv)

    suna[sunz > 0.05] = radeg * np.arctan2(sune[sunz > 0.05], sunn[sunz > 0.05])
    suna[suna < 0.0] = suna[suna < 0.0] + 360.0
    return sunz, suna
