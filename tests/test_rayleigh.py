from pathlib import Path
from rayleigh import Rayleigh, ray_press_wang
import numpy as np


def test_rayleigh():
    sensorDir = "goci"
    taur = np.array(
        [
            0.316799998,
            0.234699994,
            0.155399993,
            0.0936200023,
            0.0463000014,
            0.04098,
            0.0283000004,
            0.0155800004,
        ]
    )
    wave = np.array([412, 443, 490, 555, 660, 680, 745, 865])
    Fo = np.array(
        [
            173.914246,
            189.81572,
            197.264984,
            184.058105,
            152.560349,
            148.07373,
            128.289001,
            95.8194427,
        ]
    )

    solz = np.array([[38.5155678, 38.5155678], [38.5155678, 38.5155678]])
    senz = np.array([[36.6777306, 36.6777306], [36.6777306, 36.6777306]])
    airmass = 1.0 / np.cos(np.deg2rad(solz)) + 1.0 / np.cos(np.deg2rad(senz))
    delphi = np.array([[172.144882, 172.144882], [172.144882, 172.144882]])
    ws = np.array([[3.78513718, 3.78513718], [3.78513718, 3.78513718]])
    pr = np.array([[1018.85895, 1018.85895], [1018.85895, 1018.85895]])
    OCDATAROOT = Path("share")
    ray = Rayleigh(
        sensorDir, taur, wave, Fo, solz, senz, delphi, airmass, ws, pr, OCDATAROOT
    )
    assert np.isclose(ray.Lr[0, 0, 0], 8.31096649, rtol=1e-6)
