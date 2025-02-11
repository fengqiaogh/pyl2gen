from l12_parms import STDPR
import numpy as np
from rayleigh import Rayleigh
from whitecaps import whitecaps
from gas_trans import gaseous_transmittance
from ancillary import Ancillary
from libl1.getglint import getglint_iqu


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
    airmass = level1.airmass

    level1.t_sol = np.exp(-0.5 * pr[..., np.newaxis] / p0 * Tau_r / mu0)
    level1.t_sen = np.exp(-0.5 * pr[..., np.newaxis] / p0 * Tau_r / mu)

    # apply gaseous transmittance
    tg_sol = 1.0
    tg_sen = 1.0
    tg = 1.0

    tg_sol, tg_sen, tg = gaseous_transmittance(
        settings, ancillary, level1.sensorinfo, mu0, mu, tg_sol, tg_sen, tg
    )
    level1.tg_sol = tg_sol
    level1.tg_sen = tg_sen

    # white-cap radiances at TOA
    rhof = whitecaps(wave, ancillary.ws, ws_max=settings["wsmax"])
    level1.tLf = rhof * level1.t_sen * level1.t_sol * Fo * mu0 / np.pi

    # Rayleigh scattering
    rayleigh = Rayleigh(
        level1.sensorinfo.sensorDir,
        Tau_r,
        wave,
        Fo,
        solz,
        senz,
        raz,
        airmass,
        ancillary.ws,
        pr,
        settings["OCDATAROOT"],
    )
    level1.Lr = rayleigh.Lr

    # glint coefficients and approximate glint radiances
    # also add glint to polarization components

    glint_coef, glint_coef_q, glint_coef_u = getglint_iqu(
        senz,
        solz,
        raz,
        ancillary.ws,
    )
    level1.glint_coef = glint_coef
    level1.glint_coef_q = glint_coef_q
    level1.glint_coef_u = glint_coef_u
    level1.TLg = (
        glint_coef[..., np.newaxis]
        * np.exp(-(Tau_r + 0.1) * airmass[..., np.newaxis])
        * Fo
    )
    # 计算 L_q 和 L_u
    level1.L_q = glint_coef_q[..., np.newaxis] * level1.TLg
    level1.L_u = glint_coef_u[..., np.newaxis] * level1.TLg
