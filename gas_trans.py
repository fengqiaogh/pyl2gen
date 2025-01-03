import re
from l12_parms import O3_BIT, CO2_BIT, NO2_BIT, H2O_BIT

import numpy as np


def ozone_transmittance(oz, k_oz, csolz, csenz, tg_sol, tg_sen, tg):
    tau_oz = oz * k_oz
    tg_sol = tg_sol * np.exp(-(tau_oz / csolz))
    tg_sen = tg_sen * np.exp(-(tau_oz / csenz))
    tg = tg * np.exp(-tau_oz * (1.0 / csolz + 1.0 / csenz))
    return tg_sol, tg_sen, tg


def co2_transmittance():
    pass


def no2_transmittance():
    sec0 = 1.0 / csolz
    sec = 1.0 / csenz
    no2_tr200 = no2_frac * no2_tropo
    a_285 = k_no2 * (1.0 - 0.003 * (285.0 - 294.0))
    a_225 = k_no2 * (1.0 - 0.003 * (225.0 - 294.0))
    tau_to200 = a_285 * no2_tr200 + a_225 * no2_strat
    tg_sol = tg_sol * np.exp(-(tau_to200 * sec0))
    tg_sen = tg_sen * np.exp(-(tau_to200 * sec))
    tg = tg * np.exp(-(tau_to200 * (sec0 + sec)))
    return tg_sol, tg_sen, tg


def h2o_transmittance():
    t_h2o = a_h2o + wv * (
        b_h2o + wv * (c_h2o + wv * (d_h2o + wv * (e_h2o + wv * (f_h2o + wv * g_h2o))))
    )
    tg_sol = tg_sol * np.power(t_h2o, 1.0 / csolz)
    tg_sen = tg_sen * np.power(t_h2o, 1.0 / csenz)
    tg = tg * np.power(t_h2o, 1.0 / csenz + 1.0 / csolz)
    return tg_sol, tg_sen, tg


def gaseous_transmittance(settings: dict):
    if settings["gas_opt"] & O3_BIT:
        tg_sol, tg_sen, tg = ozone_transmittance()
    # if settings["gas_opt"] & O3_BIT:
    # tg_sol, tg_sen, tg = co2_transmittance()
    if settings["gas_opt"] & CO2_BIT:
        tg_sol, tg_sen, tg = no2_transmittance()
    if settings["gas_opt"] & H2O_BIT:
        tg_sol, tg_sen, tg = h2o_transmittance()

    return tg_sol, tg_sen, tg
