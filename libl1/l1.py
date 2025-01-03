from tkinter import N
import numpy as np

import atmocor1
from libl1.goci import GOCIL1
from oel_hdf4.filetype.filetype import FileType
from oel_hdf4.libnav.esdist import esdist as esdist_
from atmocor1 import atmocor1

FATAL_ERROR = 1


class Level1:
    def __init__(self):
        self.name = None
        self.format = None
        self.sensor_id = None
        self.subsensor_id = None

        # sensor characteristics
        self.iwave = None

        self.fsol = None

    def open(self):
        pass

    def read(self):
        match self.format:
            case FileType.FT_GOCIL1B:
                sensor_l1 = GOCIL1()
                sensor_l1.open(self.name)
                self.spatialResolution = "500 m"

            case _:
                print(f"readl1 - Unknown L1 input file format specifier: {self.format}")
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

        # Compute relative azimuth
        delphi = self.sena - 180.0 - self.sola
        delphi[delphi < -180.0] = delphi[delphi < -180.0] + 360.0
        delphi[delphi > 180.0] = delphi[delphi > 180.0] - 360.0
        self.delphi = delphi

        # Apply vicarious calibration
        # self.Lt = self.Lt * l1_input.gain + l1_input.offset

        # Scattering angle

        temp = np.sqrt(
            (1.0 - self.csenz * self.csenz) * (1.0 - self.csolz * self.csolz)
        ) * np.cos(np.deg2rad(self.delphi))
        scattang = np.acos(np.maximum(-self.csenz * self.csolz + temp, -1.0))
        self.scattang = np.rad2deg(scattang)

        self.datetime = sensor_l1.time_obj_utc

    def load(self, settings: dict):
        print(f"Loading land mask information from {settings["land"]}")
        print(f"Loading DEM information from {settings["demfile"]}")
        print(f"Loading ice mask file from {settings["icefile"]}")

        # 计算 day of year (DOY)
        day_of_year = self.datetime.timetuple().tm_yday

        # 计算毫秒 (ms)
        ms = self.datetime.second * 1000
        esdist = esdist_(self.datetime.year, day_of_year, ms)
        self.fsol = np.power(1.0 / esdist, 2)

        self.Fo = self.Fobar * self.fsol
        atmocor1(self, settings)
