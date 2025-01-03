import numpy as np


def whitecaps(wave, ws, ws_max):
    awhite = whitecap_spectral_shape(wave)
    ws_min = 6.33
    rhowc = np.where(
        (ws > ws_min) & (ws_max > ws_min),
        1.925e-5 * np.power((np.minimum(ws, ws_max) - ws_min), 3.00),
        0.0,
    )
    rhof = awhite * rhowc[..., np.newaxis]

    return rhof


def whitecap_spectral_shape(wave: np.ndarray) -> np.ndarray:
    awc_wav = np.array([412, 443, 490, 510, 555, 670, 765, 865])
    awc_tab = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 0.889225, 0.760046, 0.644950])
    return np.interp(wave, awc_wav, awc_tab, right=0.0)
