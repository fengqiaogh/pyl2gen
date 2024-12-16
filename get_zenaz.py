import numpy as np


def get_zenaz(pos, lon, lat):
    re = 6378.137
    f = 1 / 298.257
    omf2 = (1.0 - f) * (1.0 - f)

    xlat = np.deg2rad(lat)
    xlon = np.deg2rad(lon)

    npix = xlon.size
    up = np.zeros(shape=(npix, 3))
    ea = np.zeros(shape=(npix, 3))
    gv = np.zeros(shape=(npix, 3))
    up[:, 0] = np.cos(xlat) * np.cos(xlon)
    up[:, 1] = np.cos(xlat) * np.sin(xlon)
    up[:, 2] = np.sin(xlat)

    upxy = np.sqrt(up[:, 0] * up[:, 0] + up[:, 1] * up[:, 1])
    ea[:, 0] = -up[:, 1] / upxy
    ea[:, 1] = up[:, 0] / upxy
    ea[:, 2] = 0.0
    no = np.cross(up, ea)
    xlatg = np.arctan(np.tan(xlat) * omf2)

    gv[:, 0] = np.cos(xlatg) * np.cos(xlon)
    gv[:, 1] = np.cos(xlatg) * np.sin(xlon)
    gv[:, 2] = np.sin(xlatg)
    r = re * (1.0 - f) / np.sqrt(1.0 - (2.0 - f) * f * np.power(np.cos(xlatg), 2))
    r = np.tile(r, 3).reshape(3, npix)
    rh = pos - r.T * gv

    xlMatrix = np.array([ea, no, up])
    rl = np.matmul(
        np.transpose(xlMatrix, axes=(1, 0, 2)), np.transpose(rh, axes=(1, 0))
    )

    rl = np.row_stack([r[:, i] for i, r in enumerate(rl)])

    senz = np.rad2deg(
        np.arctan2(np.sqrt(rl[:, 0] * rl[:, 0] + rl[:, 1] * rl[:, 1]), rl[:, 2])
    )
    sena = np.rad2deg(np.arctan2(rl[:, 0], rl[:, 1]))

    sena[senz <= 0.05] = 0.0
    sena[sena < 0.0] = sena[sena < 0.0] + 360
    return senz, sena
