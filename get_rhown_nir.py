from sensorDefs import MAXWAVE_VIS


def get_rhown_eval(wave, nir_l, aw, bbw, chl):
    if wave[nir_l] < MAXWAVE_VIS:
        rhown_red()
    else:
        rhown_nir()


def rhown_nir(chl, aw, bbw, nir_s, nir_l):
    pass
