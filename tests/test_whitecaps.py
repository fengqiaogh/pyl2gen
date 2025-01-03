from whitecaps import whitecaps
import numpy as np


def test_whitecaps():
    wave = np.array([412, 443, 490, 555, 660, 680, 745, 865])
    ws = np.array([[3.78513718, 3.78513718], [3.78513718, 3.78513718]])
    rhof = whitecaps(wave, ws, ws_max=12)

    assert np.allclose(rhof, 0)
