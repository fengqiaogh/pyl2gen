import numpy as np
from oel_util.libgenutils.genutils import BAD_FLT
from oel_util.libgenutils.sensorDefs import MAXWAVE_VIS
from libl1.windex import Bindex, windex
from brdf import Brdf


class Chl:
    def __init__(self, wave: np.ndarray, fqfile, aw, bbw):
        self.wave = wave
        self.fqfile = fqfile
        self.aw = aw
        self.bbw = bbw
        self.chlmin = 0.001
        self.chlmax = 1000.0
        self.chlbad = BAD_FLT
        self.chl_min = 0.2
        self.chl_max = 30.0
        bindex = Bindex(wave)

        # chl hu
        ib1 = bindex.get(443)
        ib2 = bindex.get_555()
        ib3 = bindex.get(670)
        if ib3 < 0:
            ib3 = bindex.get(665)
        if ib3 < 0:
            ib3 = bindex.get(655)
        if ib3 < 0:
            ib3 = bindex.get(620)

        if ib1 < 0 or ib2 < 0 or ib3 < 0:
            print("chl_hu: incompatible sensor wavelengths for this algorithm")
            print(f"chl_hu: {ib1} {ib2} {ib3}")
            exit(1)
        else:
            self.ib1 = ib1
            self.ib2 = ib2
            self.ib3 = ib3
            print(
                f"chl_hu: using {wave[self.ib1]:7.2f} {wave[self.ib2]:7.2f} {wave[self.ib3]:7.2f}"
            )
        # chl oc3
        self.chloc3_wave = np.array([443, 489, 555])
        self.chloc3_coef = np.array([0.2515, -2.3798, 1.5823, -0.6372, -0.5692])
        oc3_ib1 = bindex.get(self.chloc3_wave[0])
        oc3_ib2 = bindex.get(self.chloc3_wave[1])
        oc3_ib3 = bindex.get(self.chloc3_wave[2])
        if oc3_ib1 < 0 or oc3_ib2 < 0 or oc3_ib3 < 0:
            print("chl_oc3: incompatible sensor wavelengths for this algorithm")
            exit(1)
        else:
            self.oc3_ib1 = oc3_ib1
            self.oc3_ib2 = oc3_ib2
            self.oc3_ib3 = oc3_ib3

        # get_rhown_nir
        ib6 = windex(670, wave)
        if np.abs(670 - wave[ib6]) > 50:
            print("can't find reasonable red band")
            print(
                f"looking for 670, found { wave[ib6]}",
            )
            exit(1)
        self.rhown_ib6 = ib6

        ib5 = windex(555, wave)
        if np.abs(555 - wave[ib5]) > 15:
            print("can't find reasonable green band")
            print(
                f"looking for 555, found { wave[ib5]}",
            )
            exit(1)
        self.rhown_ib5 = ib5

        ib2 = windex(443, wave)
        if np.abs(443 - wave[ib2]) > 5:
            print("can't find reasonable blue band")
            print(
                f"looking for 443, found { wave[ib2]}",
            )
            exit(1)
        self.rhown_ib2 = ib2
        self.brdf = Brdf(self.fqfile)

    def get(self, Rrs, method="oci"):
        match method:
            case "oci":
                return self.chl_oci(Rrs)
            case "ocx":
                return self.chl_ocx(Rrs)
            case "oc3":
                return self.chl_oc3(Rrs)
            case "hu":
                return self.chl_hu(Rrs)
            case _:
                raise ValueError("Method not supported")

    def chl_oci(self, Rrs):
        t1 = 0.25
        t2 = 0.35
        chl1 = self.chlbad
        chl2 = self.chlbad
        chl = self.chlbad
        chl1 = self.chl_hu(Rrs)
        if chl1 <= t1:
            chl = chl1
        else:
            chl2 = self.chl_ocx(Rrs)
        if chl1 >= t2:
            chl = chl2

        return chl

    def chl_ocx(self, Rrs):
        chl = self.chlbad
        chl = self.chl_oc3(Rrs)
        return chl

    def chl_oc3(self, Rrs):
        minrat = 0.201
        maxrat = 30
        ib1 = self.oc3_ib1
        ib2 = self.oc3_ib2
        ib3 = self.oc3_ib3
        Rrs1 = Rrs[ib1]
        Rrs2 = Rrs[ib2]
        Rrs3 = Rrs[ib3]
        minRrs = min(Rrs1, Rrs2)
        if Rrs3 > 0 and Rrs2 > 0 and minRrs > -0.001:
            rat = max(Rrs1, Rrs2) / Rrs3
            if rat > minrat and rat < maxrat:
                rat = np.log10(rat)
                chl = np.power(
                    10,
                    (
                        self.chloc3_coef[0]
                        + rat
                        * (
                            self.chloc3_coef[1]
                            + rat
                            * (
                                self.chloc3_coef[2]
                                + rat
                                * (self.chloc3_coef[3] + rat * self.chloc3_coef[4])
                            )
                        )
                    ),
                )
        else:
            pass
        chl = max(chl, self.chlmin)
        chl = min(chl, self.chlmax)
        return chl

    def chl_hu(self, Rrs):
        minrat = 0.21
        maxrat = 30.0
        w = np.array([443, 555, 670])
        a = np.array([-0.4287, 230.47])

        Rrs1 = Rrs[self.ib1]
        Rrs2 = Rrs[self.ib2]
        Rrs3 = Rrs[self.ib3]

        # We require that the red channel Rrs is valid (may be slightly negative
        # in clear water due to noise), and that Rrs in blue and green be positive.
        # Cases of negative blue radiance are likely high chl anyway.
        if Rrs3 > BAD_FLT + 1 and Rrs2 > 0.0 and Rrs1 > 0.0:
            # Rrs2 = conv_rrs_to_555(Rrs2)
            ci = np.minimum(
                Rrs2 - (Rrs1 + (w[1] - w[0]) / (w[2] - w[0]) * (Rrs3 - Rrs1)), 0.0
            )
            chl = np.power(10, a[0] + a[1] * ci)
            chl = max(chl, self.chlmin)
            chl = min(chl, self.chlmax)
        return chl

    def get_rhown_eval(self, wave, Rrs, nir_s, nir_l, solz, senz, phi, chl, rhown):
        if wave[nir_l] < MAXWAVE_VIS:
            return self.rhown_red()
        else:
            return self.rhown_nir(Rrs, nir_s, nir_l, solz, senz, phi, chl, rhown)

    def rhown_nir(self, Rrs, nir_s, nir_l, solz, senz, phi, chl, rhown):
        ib2 = self.rhown_ib2
        ib5 = self.rhown_ib5
        ib6 = self.rhown_ib6

        Rrs2 = Rrs[ib2]
        Rrs5 = Rrs[ib5]
        Rrs6 = Rrs[ib6]
        if Rrs6 <= 0.0:
            rhown[nir_s : nir_l + 1] = 0.0
            return rhown

        # 将 chl 限制在 chl_min 和 chl_max 之间
        chl = max(min(chl, self.chl_max), self.chl_min)

        # NOMAD fit of apg670 to chl
        apg6 = np.exp(np.log(chl) * 0.9389 - 3.7589)
        apg6 = min(max(apg6, 0.0), 0.5)

        # Compute total absorption at 670
        aw6 = self.aw[ib6]
        a6 = aw6 + apg6

        # Go below...
        Rrs2 = above_to_below(Rrs2)
        Rrs5 = above_to_below(Rrs5)
        Rrs6 = above_to_below(Rrs6)

        foq = self.brdf.foqint_morel(self.wave, solz, senz, phi, chl)

        # Compute the backscatter slope ala Lee
        if Rrs5 > 0.0 and Rrs2 > 0.0:
            eta = 2.0 * (1.0 - 1.2 * np.exp(-0.9 * (Rrs2 / Rrs5)))
            eta = min(max(eta, 0.0), 1.0)

        # Compute total backscatter at 670
        Rrs6_star = Rrs6 / foq[ib6]
        bbp6 = (Rrs6_star * a6 / (1.0 - Rrs6_star)) - self.bbw[ib6]

        # Compute normalized water-leaving reflectance at each NIR wavelength
        for ib in range(nir_s, nir_l + 1):
            if ib == ib6:
                a = a6
            else:
                a = self.aw[ib]

            # Translate bb to NIR wavelength
            bb = bbp6 * pow((self.wave[ib6] / self.wave[ib]), eta) + self.bbw[ib]

            # Remote-sensing reflectance
            salbedo = bb / (a + bb)

            Rrs_nir = foq[ib6] * salbedo

            # Normalized water-leaving reflectance
            Rrs_nir = below_to_above(Rrs_nir)
            rhown[ib] = np.pi * Rrs_nir

        return rhown


def above_to_below(Rrs):
    """
    将 Rrs[0+] 转换为 Rrs[0-]
    :param Rrs: 输入的 Rrs[0+] 值
    :return: 转换后的 Rrs[0-] 值
    """
    return Rrs / (0.52 + 1.7 * Rrs)


def below_to_above(Rrs):
    """
    将 Rrs[0-] 转换为 Rrs[0+]
    :param Rrs: 输入的 Rrs[0-] 值
    :return: 转换后的 Rrs[0+] 值
    """
    return (Rrs * 0.52) / (1 - 1.7 * Rrs)
