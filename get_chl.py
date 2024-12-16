import numpy as np


class Chl:
    chlmin = 0.00100000005
    chlmax = 1000

    def __init__(self):
        pass

    def get_chl(self, Rrs, method="default"):
        match method:
            case "x":
                return self.chl_oci(Rrs)
            case _:
                return self.chl_oci(Rrs)

    def chl_oci(self, Rrs):
        t1 = 0.25
        t2 = 0.35
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
        ib1 = 1
        ib2 = 3
        ib3 = 4
        Rrs1 = Rrs[ib1]
        Rrs2 = Rrs[ib2]
        Rrs3 = Rrs[ib3]

        ci = min(Rrs2 - (Rrs1 + (w[1] - w[0]) / (w[2] - w[0]) * (Rrs3 - Rrs1)), 0.0)
        chl = np.power(10, a[0] + a[1] * ci)

        return chl
