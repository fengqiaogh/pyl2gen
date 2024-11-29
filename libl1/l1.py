import numpy as np

from esdist import esdist
from oel_hdf4.filetype.filetype import FileType

FATAL_ERROR = 1


class L1str:
    def __init__(self):
        self.length = ""  # number of bytes allocated to data block */
        self.npix = ""

        self.iscan = ""
        self.detnum = ""
        self.mside = ""

        # scan-time-specific data */
        self.scantime = ""
        self.fsol = ""

        self.is_l2 = (
            ""  # *< Lt values are actually (above water?) reflectance skip atmocor */
        )

        # scan attributes */

        self.tilt = ""
        self.alt = ""  # altitude of sensor

        # All parameters below are scan-length dependent */

        # sensor band-pass-specific data */

        self.data = ""  # points to start of variable-length data block */

        self.nobs = ""
        self.lon = ""
        self.lat = ""
        self.solz = ""
        self.sola = ""
        self.senz = ""
        self.sena = ""
        self.Lt = ""

        self.Ltir = ""
        self.Bt = ""

        self.delphi = ""
        self.csolz = ""
        self.csenz = ""
        self.pixnum = ""
        self.slot = ""  # *< slot number                                */
        self.alpha = ""
        self.scattang = ""

        self.ws = ""
        self.wd = ""
        self.mw = ""
        self.zw = ""
        self.pr = ""
        self.oz = ""
        self.wv = ""
        self.rh = ""
        self.no2_tropo = ""
        self.no2_strat = ""
        self.no2_frac = ""
        self.sfcp = ""
        self.sfcrh = ""
        self.sfct = ""
        self.icefr = ""
        self.height = ""
        self.dem = ""

        # TODO: can get rid of this.  Only used in setanc.c
        self.ancqc = ""

        self.ssttype = ""  # per pixel - reference type or climatology */

        self.flags = ""
        self.mask = ""  # this group of params is the flags expanded into a byte
        self.hilt = ""
        self.cloud = ""
        self.glint = ""
        self.land = ""
        self.swater = ""
        self.ice = ""
        self.solzmax = ""
        self.senzmax = ""
        self.stlight = ""
        self.absaer = ""
        self.navfail = ""
        self.navwarn = ""

        self.filter = ""

        self.t_h2o = ""
        self.t_o2 = ""
        self.tg_sol = ""
        self.tg_sen = ""
        self.t_sol = ""
        self.t_sen = ""
        self.rhof = ""
        self.tLf = ""
        self.Lr = ""
        self.L_q = ""
        self.L_u = ""
        self.polcor1 = ""
        self.dpol = ""
        self.TLg = ""
        self.rhos = ""
        self.glint_coef = ""
        self.cloud_albedo = ""
        self.aerindex = ""
        self.sstref = ""
        self.sssref = ""
        self.sw_n = ""
        self.sw_a = ""
        self.sw_bb = ""
        self.sw_a_avg = ""
        self.sw_bb_avg = ""
        self.rho_cirrus = ""

        # TODO: move MERIS L1 to private_data pointer in l1rec
        # for MERIS L1 */
        self.pixdet = ""  # detector index of pixel */
        self.radcor = ""  # smile correction */

        self.Fo = ""

        # TODO: this needs to go into private_data pointer in filehandle
        # for VIIRS unaggregated and superscan */
        self.scn_fmt = ""  # scan format of data, 0 std, else unaggregated */
        self.margin_s = ""  # extra scan margin beyond actual samples */

        self.l1file = ""

        # pointer to data needed by specific readers so far just meris, hawkeye
        self.private_data = ""

        # geometry per band
        self.geom_per_band = ""

        # added ancillary data, for CHIMAERA profiles, etc
        self.anc_add = ""

        # ancillary aerosol information from MERRA-2
        self.anc_aerosol = ""

        # cloud processing data
        self.cld_dat = ""

        # uncertainty record
        self.uncertainty = ""

    def read(self, l1file, l1rec):
        match l1file.format:
            case FileType.FT_GOCIL1B:
                # readl1_goci(l1file, l1rec, 0)
                pass
            case _:
                print(
                    f"readl1 - Unknown L1 input file format specifier: {l1file.format}"
                )

    def load(self, l1_input="", l1file=""):

        yr = year
        dy = day
        ms = int(sec * 1.0e3)

        # Apply vicarious calibration
        l1rec.Lt = l1rec.Lt * l1_input.gain + l1_input.offset

        # Compute relative azimuth
        delphi = self.sena - 180.0 - self.sola
        delphi[delphi < -180.0] = delphi[delphi < -180.0] + 360.0
        delphi[delphi > 180.0] = delphi[delphi > 180.0] - 360.0
        self.delphi = delphi

        # Precompute frequently used trig relations
        self.csolz = np.cos(np.deg2rad(self.solz))
        self.csenz = np.cos(np.deg2rad(self.senz))

        # Scattering angle
        """
        temp = sqrt((1.0 - l1rec->csenz[ip] * l1rec->csenz[ip])*(1.0 - l1rec->csolz[ip] * l1rec->csolz[ip]))
                * cos(l1rec->delphi[ip] / radeg);
        l1rec->scattang[ip] = acos(MAX(-l1rec->csenz[ip] * l1rec->csolz[ip] + temp, -1.0)) * radeg;

        """
