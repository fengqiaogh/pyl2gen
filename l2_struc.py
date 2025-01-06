from libl1.l1 import Level1


class Level2:
    def __init__(self):
        self.l1rec = Level1()
        self.length = ""
        self.data = ""

        # var[npix]
        self.num_iter = ""
        self.aermodmin = ""
        self.aermodmax = ""
        self.aermodmin2 = ""
        self.aermodmax2 = ""

        self.chl = ""
        self.eps = ""  # NIR aerosol reflectance ratio (single scattering)
        self.aerratio = ""
        self.aerratio2 = ""
        self.aerindex = ""

        # var[npixnbands]
        self.taua = ""  # aerosol optical thickness
        self.La = ""  # aerosol radiance1
        self.Lw = ""  # water-leaving radiance
        self.nLw = ""  # normalized water-leaving radiance
        self.nLw_unc = ""
        self.brdf = ""  # bi-direction reflectance function
        self.Rrs = ""  # Remote sensing reflectance
        self.Rrs_unc = ""
        self.chi2 = ""  # chi square from the spectral match in mbac
        self.chl_unc = ""
        self.outband_correction = ""  # square bandpass correction for Rrs
        self.a = ""  # absoprtion coefficient
        self.bb = ""  # backscattering coefficient

        # allocated or set later
        self.bindx = ""
        self.sst = ""
        self.Rrs_raman = ""
        self.tgrec = ""
