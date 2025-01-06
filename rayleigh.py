import numpy as np
from scipy.interpolate import RegularGridInterpolator
import h5py

NSOL = 45
NSEN = 41
NORDER = 3
NWIND = 8
STDPR = 1013.25


class Rayleigh:
    def __init__(
        self,
        sensor,
        taur,
        wave: np.ndarray,
        Fo: np.ndarray,
        solz: np.ndarray,
        senz: np.ndarray,
        delphi: np.ndarray,
        airmass: np.ndarray,
        ws: np.ndarray,
        pr: np.ndarray,
        OCDATAROOT,
    ):

        self.Lr = np.zeros((ws.shape[0], ws.shape[1], wave.size))
        self.L_q = np.zeros((ws.shape[0], ws.shape[1], wave.size))
        self.L_u = np.zeros((ws.shape[0], ws.shape[1], wave.size))

        sigma_m = 0.0731 * np.sqrt(ws)

        fac = ray_press_wang(taur, airmass, pr)

        for idx, wa in enumerate(wave):

            file = OCDATAROOT.joinpath(
                f"{sensor}", "rayleigh", f"rayleigh_{sensor}_{wa}_iqu.h5"
            )

            print(f"Loading Rayleigh LUT {file}")
            with h5py.File(file) as f:
                senz_tab = f["senz"][()]
                solz_tab = f["solz"][()]
                sigma_tab = f["sigma"][()]

                # I component
                i_ray_tab = f["i_ray"][()]

            tmp = np.zeros((NWIND, NSOL, NSEN))
            ray_i = np.sum(
                [
                    RegularGridInterpolator(
                        (sigma_tab, solz_tab, senz_tab), i_ray_tab[:, :, m, :]
                    )(np.stack([sigma_m, solz, senz], axis=-1))
                    * np.cos(np.deg2rad(delphi) * m)
                    for m in range(NORDER)
                ],
                axis=0,
            )

            self.Lr[:, :, idx] = Fo[idx] * ray_i * fac[:, :, idx]


def ray_press_wang(taur, airmass, pr) -> np.ndarray:
    p0 = STDPR
    x = (
        (
            -(0.6543 - 1.608 * taur)
            + (0.8192 - 1.2541 * taur) * np.log(airmass[..., np.newaxis])
        )
        * taur
        * airmass[..., np.newaxis]
    )
    return (1.0 - np.exp(-x * pr[..., np.newaxis] / p0)) / (1.0 - np.exp(-x))
