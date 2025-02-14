import numpy as np
from ancillary import Ancillary
import atmocor1
from libl1.goci import GOCIL1
from oel_hdf4.filetype.filetype import FileType
from oel_hdf4.libnav.esdist import esdist as esdist_
from atmocor1 import atmocor1
from polcor import polcor
from libl1.setflags import Flag
from libl1.windex import Bindex

FATAL_ERROR = 1
OFF = 0
ON = 1
BANDW = 10


class Level1:
    def __init__(self):
        self.name = None
        self.format = None
        self.wave = None

        # sensor information
        self.sensorinfo = None

        self.fsol = None
        self.Fo = None

        self.Lt = None

        # Rayleigh
        self.Lr = None

        # white-cap radiances at TOA
        self.tLf = None

        # add surface reflectance
        self.rhos = None
        self.airmass = None

        # solar angle
        self.solz = None
        self.csolz = None
        self.sola = None

        # satellite angle
        self.senz = None
        self.csenz = None
        self.sena = None

        # land mask
        self.land = None

        # 气体透过率
        self.tg_sol = None
        self.tg_sen = None
        self.t_sol = None
        self.t_sen = None
        self.t_o2 = 1.0
        self.t_h2o = 1.0

        self.glint_coef = None
        self.glint_coef_q = None
        self.glint_coef_u = None

        self.bindex = Bindex()

    def read(self, sline, eline, spixl, epixl):

        match self.format:
            case FileType.FT_GOCIL1B:
                sensor_l1 = GOCIL1()
                sensor_l1.open(self.name, sline, eline, spixl, epixl)
                self.spatialResolution = "500 m"
                self._process_sensor_data(sensor_l1)

            case _:
                print(f"readl1 - Unknown L1 input file format specifier: {self.format}")

    def _process_sensor_data(self, sensor_l1):

        self.npix = sensor_l1.npixels
        self.nscan = sensor_l1.nscans
        self.bands = sensor_l1.nbands
        self.lat = sensor_l1.latitudes
        self.lon = sensor_l1.longitudes
        self.solz = sensor_l1.solz
        self.csolz = np.cos(np.deg2rad(self.solz))
        self.sola = sensor_l1.sola
        self.senz = sensor_l1.senz
        self.csenz = np.cos(np.deg2rad(self.senz))
        self.sena = sensor_l1.sena
        self.airmass = 1.0 / self.csolz + 1.0 / self.csenz

        self.Lt = sensor_l1.Lt
        self._compute_relative_azimuth()
        self._compute_scattering_angle()
        self.datetime = sensor_l1.time_obj_utc

        self.bindex.set(self.wave, BANDW)

    def _compute_relative_azimuth(self):
        # Compute relative azimuth
        delphi = self.sena - 180.0 - self.sola
        delphi[delphi < -180.0] = delphi[delphi < -180.0] + 360.0
        delphi[delphi > 180.0] = delphi[delphi > 180.0] - 360.0
        self.delphi = delphi

    def _compute_scattering_angle(self):
        # Scattering angle
        temp = np.sqrt(
            (1.0 - self.csenz * self.csenz) * (1.0 - self.csolz * self.csolz)
        ) * np.cos(np.deg2rad(self.delphi))
        scattang = np.acos(np.maximum(-self.csenz * self.csolz + temp, -1.0))
        self.scattang = np.rad2deg(scattang)

    def load(self, settings: dict, ancillary: Ancillary):
        # 计算 day of year (DOY)
        day_of_year = self.datetime.timetuple().tm_yday

        # 计算毫秒 (ms)
        ms = self.datetime.second * 1000
        esdist = esdist_(self.datetime.year, day_of_year, ms)
        self.fsol = np.power(1.0 / esdist, 2)

        self.Fo = self.sensorinfo.F0 * self.fsol

        # Apply vicarious calibration
        self.Lt = self.Lt * settings["gain"] + settings["offset"]

        atmocor1(self, settings, ancillary)

        # set polarization correction
        polcor()

        # add surface reflectance
        self.get_rhos()

        # set masks and flags
        flag = Flag(self.sensorinfo.Lambda, settings, self.lon, self.lat)
        flag.set(self)

    def get_rhos(self):
        mu0 = self.csolz[..., np.newaxis]
        self.rhos = (
            np.pi
            / self.Fo
            / mu0
            * (self.Lt / self.tg_sol / self.tg_sen - self.Lr)
            / self.t_sol
            / self.t_sen
            / self.t_o2
            / self.t_h2o
        )
