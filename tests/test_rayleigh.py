from rayleigh import Rayleigh, ray_press_wang
import numpy as np


def test_ray_press_wang():
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
    airmass = np.array([[2.52492547, 2.52492547], [2.52492547, 2.52492547]])
    pr = np.array([[1018.85895, 1018.85895], [1018.85895, 1018.85895]])
    fac = ray_press_wang(taur, airmass, pr)
    assert np.allclose(
        fac[0, 0, :],
        np.array(
            [
                1.00500631,
                1.00519776,
                1.00534809,
                1.00543988,
                1.00549495,
                1.00550032,
                1.00551236,
                1.00552344,
            ]
        ),
        rtol=1e-6,
    )
