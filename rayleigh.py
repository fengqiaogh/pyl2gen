import os
import sys

import numpy as np
from pyhdf import SD
from scipy.interpolate import RegularGridInterpolator

from libpiutils.sensorInfo import SensorInfo

NSOL = 45
NSEN = 41
NORDER = 3
NWIND = 8
STDPR = 1013.25


class Rayleigh:
    def __init__(self, l1rec, pol_opt=0):

        self.sensorID = l1rec.l1file.sensorID
        self.nwave = l1rec.l1file.nbands
        self.taur = l1rec.l1file.Tau_r
        self.Fo = l1rec.Fo
        self.pr = l1rec.pr
        self.ws = l1rec.ws
        self.r_solz = l1rec.solz
        self.r_senz = l1rec.senz
        self.r_phi = l1rec.delphi
        self.r_csolz = l1rec.csolz
        self.r_csenz = l1rec.csenz
        npix = len(self.r_solz)

        OCDATAROOT = os.environ.get("OCDATAROOT", None)
        if OCDATAROOT == None:
            print("OCDATAROOT environment variable is not defined.")
            sys.exit(1)

        sigma_m = 0.0731 * np.sqrt(self.ws)
        airmass = 1.0 / self.r_csolz + 1.0 / self.r_csenz
        global sen_info
        sensorstr = sen_info.sensorName.lower()
        wave = sen_info.Lambda
        for iw in range(0, self.nwave):
            fac = self.ray_press_wang(self.taur[iw], airmass, self.pr)
            wavestr = str(wave[iw])
            file = os.path.join(
                OCDATAROOT,
                f"{sen_info.sensorDir}",
                "rayleigh",
                f"rayleigh_{sensorstr}_{wavestr}_iqu.hdf",
            )

            print(f"Loading Rayleigh LUT {file}")
            f = SD.SD(file, SD.SDC.READ)

            ray_sen = f.select("senz")[:]
            ray_sol = f.select("solz")[:]
            ray_sigma = f.select("sigma")[:]

            # I component
            ray_for_i = f.select("i_ray")[:]

            tmp = np.zeros((NWIND, NSOL, NSEN))
            ray_i = np.zeros(npix)

            for m in range(0, NORDER):
                tmp[:, :, :] = ray_for_i[:, :, m, :]
                interpolater = RegularGridInterpolator(
                    (ray_sigma, ray_sol, ray_sen), tmp
                )
                target_points = np.array([sigma_m, self.r_solz, self.r_senz])
                interpolated_values = interpolater(target_points.T)
                ray_i = ray_i + interpolated_values * np.cos(np.deg2rad(self.r_phi) * m)
            l1rec.Lr[:, iw] = self.Fo[iw] * ray_i * fac

            if pol_opt > 0:
                # Q component
                self.ray_for_q = f.select("q_ray")[:]
                # U component
                self.ray_for_u = f.select("u_ray")[:]

                for m in range(0, NORDER):
                    ray_q = ""

                for m in range(0, NORDER):
                    ray_u = ""

                l1rec.L_q[:, iw] = self.Fo[iw] * ray_q * fac
                l1rec.L_u[:, iw] = self.Fo[iw] * ray_u * fac
            f.end()

    def ray_press_wang(self, taur, airmass, pr):
        p0 = STDPR
        x = (
            (-(0.6543 - 1.608 * taur) + (0.8192 - 1.2541 * taur) * np.log(airmass))
            * taur
            * airmass
        )
        return (1.0 - np.exp(-x * pr / p0)) / (1.0 - np.exp(-x))
