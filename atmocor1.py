from l12_parms import STDPR
import numpy as np
from rayleigh import Rayleigh
from whitecaps import whitecaps
from gas_trans import gaseous_transmittance


def atmocor1(level1, settings: dict):
    p0 = STDPR
    sensorID = level1.sensorID
    nwave = level1.nbands

    solz = level1.solz
    senz = level1.senz
    raz = level1.delphi
    mu0 = level1.csolz
    mu = level1.csenz

    ws = level1.ws
    pr = level1.pr
    wv = level1.wv

    Fo = level1.Fo
    Tau_r = level1.Tau_r

    airmass = 1.0 / mu0 + 1.0 / mu

    # apply gaseous transmittance
    tg_sol, tg_sen, tg = gaseous_transmittance()

    # white-cap radiances at TOA
    rhof = whitecaps(wave, ws, ws_max=settings["wsmax"])

    tLf = rhof * t_sen * t_sol * Fo * mu0 / np.pi

    # Rayleigh scattering
    rayleigh = Rayleigh()
