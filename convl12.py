from atmocor2 import atmocor2
from aerosol import Aerosol
import numpy as np
from chl import Chl


def convl12(level1, level2, ancillary, settings):
    sensorID = level1.sensorinfo.sensorId
    wave = level1.sensorinfo.Lambda
    OCDATAROOT = settings["OCDATAROOT"]
    aermodfile = settings.get("aermodfile", "").replace(
        "$OCDATAROOT", OCDATAROOT.as_posix()
    )
    aermodels = settings["aermodels"]
    naermodels = settings["naermodels"]
    aer_opt = settings["aer_opt"]

    aer_iter_max = settings["aer_iter_max"]
    nir_s = np.argwhere(wave == settings["aer_wave_short"])[0][0]
    nir_l = np.argwhere(wave == settings["aer_wave_long"])[0][0]
    aer_base = np.argwhere(wave == settings["aer_wave_base"])[0][0]
    rhoamin = settings["rhoamin"]

    if nir_s < 0 or nir_l < 0:
        raise ValueError(
            f"Aerosol selection bands {settings["aer_wave_short"]} and {settings["aer_wave_long"]} not available for this sensor."
        )
    if nir_l < nir_s:
        raise ValueError(
            f"Invalid aerosol selection bands: long ({settings["aer_wave_long"]}) must be greater than short ({settings["aer_wave_short"]})."
        )

    if wave[nir_s] < 600:
        raise ValueError(f"Aerosol selection band(s) must be greater than 600nm.")

    aerosol = Aerosol(
        sensorID,
        wave,
        aermodfile,
        aermodels,
        naermodels,
        aer_opt,
        nir_s,
        nir_l,
        aer_iter_max,
        rhoamin,
    )

    print(f"Aerosol selection bands {wave[nir_s]} and {wave[nir_l]}")

    match aer_opt:
        case -17:  # AERRHMSEPS
            want_nirLw = 1
            print("NIR correction enabled --> for multi-scattering epsilon.")

        case _:
            if settings["aer_rrs_short"] >= 0.0 and settings["aer_rrs_long"] >= 0.0:
                want_nirRrs = 1
                aer_iter_min = 3
                aer_iter_max = settings["aer_iter_max"]
                print("NIR correction via input Rrs enabled.")

    if want_nirLw or want_nirRrs:
        red = np.argwhere(wave == 670)
        if red.size == 0:
            red = np.argwhere(wave == 680)
            if red.size == 0:
                red = np.argwhere(wave == 620)
                if red.size == 0:
                    red = np.argwhere(wave == 765)
                    if red.size == 0:
                        red = np.argwhere(wave == 655)
                        if red.size == 0:
                            red = np.argwhere(wave == 664)  # added for MSI
                            if red.size == 0:
                                print("can't find red band")
                                exit(1)

        green = np.argwhere(wave == 550)
        if green.size == 0:
            green = np.argwhere(wave == 555)
            if green.size == 0:
                green = np.argwhere(wave == 560)
                if green.size == 0:
                    green = np.argwhere(wave == 565)
                    if green.size == 0:
                        print("can't find green band")
                        exit(1)

    bindex = level1.bindex
    chlorophyll = Chl(wave, bindex)
    # 在这里开始逐像素运行，方便以后并行处理
    # t_sol = level1.t_sol
    # t_sen = level1.t_sen
    # tg = level1.tg
    # Fo = level1.Fo
    Tau_r = level1.sensorinfo.Tau_r
    Fo = level1.Fo

    # Test
    solz = 39.2826195
    senz = 37.2348137
    phi = 172.361465

    t_sol = np.array(
        [
            0.814087629,
            0.858659685,
            0.904027164,
            0.94102627,
            0.970386326,
            0.973743975,
            0.981793582,
            0.989935458,
        ]
    )
    t_sen = np.array(
        [
            0.818755329,
            0.862304389,
            0.906566083,
            0.942617536,
            0.971197546,
            0.974464357,
            0.982295096,
            0.990213811,
        ]
    )
    tg = np.array(
        [
            0.99287581716725404,
            0.99223465514819775,
            0.98180457602387106,
            0.93426234230524396,
            0.96258138631671808,
            0.97531594313150993,
            0.99266292802566503,
            0.99848538440915779,
        ]
    )
    tg_sol = np.array(
        [
            0.9963813392104125,
            0.9960550402329148,
            0.99073231105634008,
            0.9661094215443865,
            0.9808486851966226,
            0.98740692655382689,
            0.99627300767642979,
            0.99923173188838366,
        ]
    )
    pr = 1018.42242
    Lt = np.array(
        [
            9.47686958,
            8.59582138,
            7.21518755,
            5.10522747,
            2.12941766,
            1.86647546,
            1.0473057,
            0.567476273,
        ]
    )
    Lr = np.array(
        [
            8.34885883,
            6.92787933,
            4.85848713,
            2.74680662,
            1.12013352,
            0.960908532,
            0.572968483,
            0.234133989,
        ]
    )
    tLf = np.zeros(8)

    rh = 70.1901093
    fsol = 1.0039427256657802
    Fobar = np.array(
        [
            173.231003,
            189.070007,
            196.490005,
            183.335007,
            151.960999,
            147.492004,
            127.785004,
            95.4430008,
        ]
    )
    num_iter, eps, Rrs, chl = atmocor2(
        Lt,
        Lr,
        tLf,
        solz,
        senz,
        phi,
        Tau_r,
        t_sol,
        t_sen,
        tg,
        tg_sol,
        Fo,
        rh,
        pr,
        aerosol,
        fsol,
        Fobar,
        chlorophyll,
    )
