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
        csolz: np.ndarray,
        csenz: np.ndarray,
        delphi: np.ndarray,
        ws: np.ndarray,
        pr: np.ndarray,
        settings: dict,
    ):

        self.Lr = np.zeros((ws.size[0], ws.size[1], wave.size))
        self.L_q = np.zeros((ws.size[0], ws.size[1], wave.size))
        self.L_u = np.zeros((ws.size[0], ws.size[1], wave.size))

        sigma_m = 0.0731 * np.sqrt(ws)
        airmass = 1.0 / csolz + 1.0 / csenz
        fac = ray_press_wang(taur, airmass, pr)
        for idx, wa in enumerate(wave):

            OCDATAROOT = settings["OCDATAROOT"]
            file = OCDATAROOT.joinpath(
                f"{sensor}", "rayleigh", f"rayleigh_{sensor}_{wa}_iqu.h5"
            )

            print(f"Loading Rayleigh LUT {file}")
            with h5py.File(file) as f:
                senz_tab = f.select("senz")[()]
                solz_tab = f.select("solz")[()]
                sigma_tab = f.select("sigma")[()]

                # I component
                i_ray_tab = f.select("i_ray")[()]

            tmp = np.zeros((NWIND, NSOL, NSEN))

            for m in range(0, NORDER):
                tmp[:, :, :] = i_ray_tab[:, :, m, :]
                interpolater = RegularGridInterpolator(
                    (sigma_tab, solz_tab, senz_tab), tmp
                )
                target_points = np.array([sigma_m, csolz, csenz])
                interpolated_values = interpolater(target_points.T)
                ray_i = ray_i + interpolated_values * np.cos(np.deg2rad(delphi) * m)
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
