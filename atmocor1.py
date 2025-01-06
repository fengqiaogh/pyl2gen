from l12_parms import STDPR
import numpy as np
from rayleigh import Rayleigh
from whitecaps import whitecaps
from gas_trans import gaseous_transmittance
from ancillary import Ancillary
from libl1.setflags import Flag


def atmocor1(level1, settings: dict, ancillary: Ancillary):
    p0 = STDPR
    pr = ancillary.pr
    Tau_r = level1.sensorinfo.Tau_r
    wave = level1.sensorinfo.Lambda
    Fo = level1.Fo
    solz = level1.solz
    senz = level1.senz
    raz = level1.delphi
    mu0 = level1.csolz[..., np.newaxis]
    mu = level1.csenz[..., np.newaxis]

    t_sol = np.exp(-0.5 * pr[..., np.newaxis] / p0 * Tau_r / mu0)
    t_sen = np.exp(-0.5 * pr[..., np.newaxis] / p0 * Tau_r / mu)

    # apply gaseous transmittance
    tg_sol = 1.0
    tg_sen = 1.0
    tg = 1.0
    t_h2o = 1.0
    t_o2 = 1.0
    tg_sol, tg_sen, tg = gaseous_transmittance(
        settings, ancillary, level1.sensorinfo, mu0, mu, tg_sol, tg_sen, tg
    )

    # white-cap radiances at TOA
    rhof = whitecaps(wave, ancillary.ws, ws_max=settings["wsmax"])
    level1.tLf = rhof * t_sen * t_sol * Fo * mu0 / np.pi

    # Rayleigh scattering
    rayleigh = Rayleigh(
        level1.sensorinfo.sensorDir,
        Tau_r,
        wave,
        Fo,
        solz,
        senz,
        raz,
        level1.airmass,
        ancillary.ws,
        pr,
        settings["OCDATAROOT"],
    )
    level1.Lr = rayleigh.Lr

    # glint coefficients and approximate glint radiances
    # also add glint to polarization components
    """
    glint_coef, glint_coef_q, glint_coef_u = getglint_iqu(
        senz,
        solz,
        raz,
        ancillary.ws,
    )
    """

    # add surface reflectance
    # get_rhos()
    level1.rhos = (
        np.pi
        / Fo
        / mu0
        * (level1.Lt / tg_sol / tg_sen - level1.Lr)
        / t_sol
        / t_sen
        / t_o2
        / t_h2o
    )

    # set masks and flags
    flag = Flag()
    flag.set()
