import numpy as np


def get_zenaz(pos, lon, lat):
    re = 6378.137
    f = 1 / 298.257
    omf2 = (1.0 - f) * (1.0 - f)

    xlat = np.deg2rad(lat)
    xlon = np.deg2rad(lon)

    nrow, ncol = lat.shape
    up = np.zeros(shape=(nrow, ncol, 3))
    ea = np.zeros(shape=(nrow, ncol, 3))
    gv = np.zeros(shape=(nrow, ncol, 3))
    up[:, :, 0] = np.cos(xlat) * np.cos(xlon)
    up[:, :, 1] = np.cos(xlat) * np.sin(xlon)
    up[:, :, 2] = np.sin(xlat)

    upxy = np.sqrt(up[:, :, 0] * up[:, :, 0] + up[:, :, 1] * up[:, :, 1])
    ea[:, :, 0] = -up[:, :, 1] / upxy
    ea[:, :, 1] = up[:, :, 0] / upxy
    ea[:, :, 2] = 0.0
    no = np.cross(up, ea)
    xlatg = np.arctan(np.tan(xlat) * omf2)

    gv[:, :, 0] = np.cos(xlatg) * np.cos(xlon)
    gv[:, :, 1] = np.cos(xlatg) * np.sin(xlon)
    gv[:, :, 2] = np.sin(xlatg)
    r = re * (1.0 - f) / np.sqrt(1.0 - (2.0 - f) * f * np.power(np.cos(xlatg), 2))
    r = np.tile(r, 3)
    r = np.reshape(r, (nrow, ncol, 3))
    rh = pos - r * gv

    xlMatrix = np.array([ea, no, up])
    xlMatrix = np.transpose(xlMatrix, axes=(1, 2, 0, 3))  # shape (nrow, ncol, 3, 3)

    # Reshape rh to (nrow, ncol, 3, 1)
    rh = np.expand_dims(rh, axis=-1)

    # Perform matrix multiplication
    rl = np.matmul(xlMatrix, rh)  # shape (nrow, ncol, 3, 3)

    # Remove the last dimension to get the result in shape (nrow, ncol, 3)
    rl = np.squeeze(rl, axis=-1)  # shape (nrow, ncol, 3)

    senz = np.rad2deg(
        np.arctan2(
            np.sqrt(rl[:, :, 0] * rl[:, :, 0] + rl[:, :, 1] * rl[:, :, 1]), rl[:, :, 2]
        )
    )
    sena = np.rad2deg(np.arctan2(rl[:, :, 0], rl[:, :, 1]))

    sena[senz <= 0.05] = 0.0
    sena[sena < 0.0] = sena[sena < 0.0] + 360
    return senz, sena
