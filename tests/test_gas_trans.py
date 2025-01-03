from gas_trans import ozone_transmittance
import numpy as np


def test_ozone_transmittance():
    oz = 0.268703669
    k_oz = np.array(
        [
            0.000472899992,
            0.002997,
            0.0220400002,
            0.0971999988,
            0.0552100018,
            0.036150001,
            0.0106600001,
            0.00219799997,
        ]
    )
    csolz = 0.782438993
    csenz = 0.802007854
    tg_sol, tg_sen, tg = ozone_transmittance(
        oz, k_oz, csolz, csenz, tg_sol=1, tg_sen=1, tg=1
    )
    assert np.isclose(tg_sol[0], 0.999837637)
    assert np.isclose(tg_sen[0], 0.999841571)
    assert np.isclose(tg[0], 0.999679208)
