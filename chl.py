import numpy as np
from oel_util.libgenutils.genutils import BAD_FLT


class Chl:
    def __init__(self, wave: np.ndarray, bindex):
        self.chlmin = 0.001
        self.chlmax = 1000.0
        self.chlbad = BAD_FLT
        self.minrat = 0.21
        self.maxrat = 30.0

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
        chl = self.chl_oc3(Rrs)
        return chl

    def chl_oc3(self, Rrs):
        a = np.array([0.25150001, -2.37980008, 1.58229995, -0.637199998, -0.569199979])
        minrat = 0.201
        maxrat = 30
        ib1 = 1
        ib2 = 2
        ib3 = 3
        Rrs1 = Rrs[ib1]
        Rrs2 = Rrs[ib2]
        Rrs3 = Rrs[ib3]
        minRrs = min(Rrs1, Rrs2)
        if Rrs3 > 0 and Rrs2 > 0 and minRrs > -0.001:
            rat = max(Rrs1, Rrs2) / Rrs3
            if rat > minrat and rat < maxrat:
                rat = np.log10(rat)
                chl = np.power(
                    10, (a[0] + rat * (a[1] + rat * (a[2] + rat * (a[3] + rat * a[4]))))
                )
        else:
            pass
        return chl

    def chl_hu(self, Rrs):
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
