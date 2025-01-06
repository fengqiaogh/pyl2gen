import numpy as np
from l12_parms import O3_BIT, CO2_BIT, NO2_BIT, H2O_BIT

from ancillary import Ancillary


def ozone_transmittance(oz, k_oz, csolz, csenz, tg_sol, tg_sen, tg):
    tau_oz = oz * k_oz
    tg_sol = tg_sol * np.exp(-(tau_oz / csolz))
    tg_sen = tg_sen * np.exp(-(tau_oz / csenz))
    tg = tg * np.exp(-tau_oz * (1.0 / csolz + 1.0 / csenz))
    return tg_sol, tg_sen, tg


def co2_transmittance():
    pass


def no2_transmittance(
    no2_strat, no2_tropo, no2_frac, k_no2, csolz, csenz, tg_sol, tg_sen, tg
):
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


def h2o_transmittance(
    wv,
    a_h2o,
    b_h2o,
    c_h2o,
    d_h2o,
    e_h2o,
    f_h2o,
    g_h2o,
    csolz,
    csenz,
    tg_sol,
    tg_sen,
    tg,
):
    t_h2o = a_h2o + wv * (
        b_h2o + wv * (c_h2o + wv * (d_h2o + wv * (e_h2o + wv * (f_h2o + wv * g_h2o))))
    )
    tg_sol = tg_sol * np.power(t_h2o, 1.0 / csolz)
    tg_sen = tg_sen * np.power(t_h2o, 1.0 / csenz)
    tg = tg * np.power(t_h2o, 1.0 / csenz + 1.0 / csolz)
    return tg_sol, tg_sen, tg


def gaseous_transmittance(
    settings: dict, ancillary: Ancillary, sensorinfo, csolz, csenz, tg_sol, tg_sen, tg
):

    if settings["gas_opt"] & O3_BIT:
        oz = ancillary.oz[..., np.newaxis]
        tg_sol, tg_sen, tg = ozone_transmittance(
            oz, sensorinfo.k_oz, csolz, csenz, tg_sol, tg_sen, tg
        )
    # if settings["gas_opt"] & O3_BIT:
    # tg_sol, tg_sen, tg = co2_transmittance()
    if settings["gas_opt"] & CO2_BIT:
        no2_strat = ancillary.no2_strat[..., np.newaxis]
        no2_tropo = ancillary.no2_tropo[..., np.newaxis]
        no2_frac = ancillary.no2_frac[..., np.newaxis]
        tg_sol, tg_sen, tg = no2_transmittance(
            no2_strat,
            no2_tropo,
            no2_frac,
            sensorinfo.k_no2,
            csolz,
            csenz,
            tg_sol,
            tg_sen,
            tg,
        )
    if settings["gas_opt"] & H2O_BIT:
        wv = ancillary.wv[..., np.newaxis]
        tg_sol, tg_sen, tg = h2o_transmittance(
            wv,
            sensorinfo.a_h2o,
            sensorinfo.b_h2o,
            sensorinfo.c_h2o,
            sensorinfo.d_h2o,
            sensorinfo.e_h2o,
            sensorinfo.f_h2o,
            sensorinfo.g_h2o,
            csolz,
            csenz,
            tg_sol,
            tg_sen,
            tg,
        )

    return tg_sol, tg_sen, tg
