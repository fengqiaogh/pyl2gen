from datetime import datetime

import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI
from l12_parms import NO2_BIT


class Ancillary:
    def __init__(self):
        self.rh = None
        self.ws = None
        self.wd = None
        self.pr = None
        self.wv = None  # precipitable water (water vapor)
        self.oz = None  # ozone

        # no2 and fraction
        self.no2_tropo = None
        self.no2_strat = None
        self.no2_frac = None

    def set(self, lat, lon, dt: datetime, settings: dict):
        OCDATAROOT = settings["OCDATAROOT"]
        self.met1 = settings.get("met1", "").replace(
            "$OCDATAROOT", OCDATAROOT.as_posix()
        )
        self.met2 = settings.get("met2", "").replace(
            "$OCDATAROOT", OCDATAROOT.as_posix()
        )
        self.met3 = settings.get("met3", "").replace(
            "$OCDATAROOT", OCDATAROOT.as_posix()
        )
        self.ozone1 = settings.get("ozone1", "").replace(
            "$OCDATAROOT", OCDATAROOT.as_posix()
        )
        self.ozone2 = settings.get("ozone2", "").replace(
            "$OCDATAROOT", OCDATAROOT.as_posix()
        )
        self.ozone3 = settings.get("ozone3", "").replace(
            "$OCDATAROOT", OCDATAROOT.as_posix()
        )
        self.no2file = settings.get("no2file", "").replace(
            "$OCDATAROOT", OCDATAROOT.as_posix()
        )

        print("Opening meteorological files.")
        print(f"  met1   = {self.met1}")
        print(f"  met2   = {self.met2}")
        print(f"  met3   = {self.met3}")
        print(f"  ozone1 = {self.ozone1}")
        print(f"  ozone2 = {self.ozone2}")
        print(f"  ozone3 = {self.ozone3}")
        print(f"  no2    = {self.no2file}")

        # relative humidity
        parmID = 5
        self.rh = get_ancillary(lat, lon, dt, self.met1, self.met2, self.met3, parmID)

        # wind speed
        parmID = 0
        zw = get_ancillary(lat, lon, dt, self.met1, self.met2, self.met3, parmID)

        parmID = 1
        mw = get_ancillary(lat, lon, dt, self.met1, self.met2, self.met3, parmID)

        # create wind speed, direction from u, v components
        u = zw
        v = mw
        ws_2 = u * u + v * v
        self.ws = np.sqrt(ws_2)
        self.wd = np.rad2deg(np.arctan2(-u, -v))

        # surface pressure
        parmID = 2
        pr = get_ancillary(lat, lon, dt, self.met1, self.met2, self.met3, parmID)
        pr[pr <= 0] = 1013.25
        pr = np.nan_to_num(pr, nan=1013.25)
        pr[pr < 900] = 900
        pr[pr > 1100] = 1100
        self.pr = pr

        # precipitable water (water vapor)
        parmID = 3
        wv = get_ancillary(lat, lon, dt, self.met1, self.met2, self.met3, parmID)

        # convert from kg/m^2 to g/cm^2
        self.wv = wv / 10

        # ozone
        jday = dt.timetuple().tm_yday
        oz = ozone_climatology(self.ozone1, jday, lon, lat)

        # convert from Dobson units to atm-cm
        self.oz = oz / 1000.0

        # no2 and fraction
        if settings["gas_opt"] & NO2_BIT:
            no2_tropo, no2_strat = no2conc(self.no2file, lon, lat, dt.month)
            self.no2_tropo = no2_tropo * 1e15
            self.no2_strat = no2_strat * 1e15

            no2_frac_file = OCDATAROOT.joinpath("common/trop_f_no2_200m.h5")

            self.no2_frac = no2_frac_(no2_frac_file, lon, lat)

    def get(self):
        pass


def no2_frac_(no2_frac_file, lon, lat):
    print(f"Opening NO2 frac file {no2_frac_file}")
    tab_lon = np.arange(-180, 180, 2)
    tab_lat = -np.arange(-90, 90, 2)
    with h5py.File(no2_frac_file) as f:
        data = f["Geophysical Data"]["f_no2_200m"][()]
    func = RGI((tab_lat, tab_lon), data)

    return func(np.stack([lat, lon], axis=-1))


def no2conc(no2file, lon, lat, month):
    print(f"Opening NO2 file {no2file}")
    tab_lon = np.arange(-180, 180, 0.25)
    tab_lat = -np.arange(-90, 90, 0.25)
    with h5py.File(no2file) as f:
        tot = f["Geophysical Data"][f"tot_no2_{month:02d}"][()]
        trop = f["Geophysical Data"][f"trop_no2_{month:02d}"][()]
    tot_func = RGI((tab_lat, tab_lon), tot)
    total = tot_func(np.stack([lat, lon], axis=-1))
    trop_func = RGI((tab_lat, tab_lon), trop)
    tropo = trop_func(np.stack([lat, lon], axis=-1))
    no2_strat = np.maximum(total - tropo, 0.0)
    no2_tropo = np.maximum(tropo, 0.0)
    return no2_tropo, no2_strat


def ozone_climatology(file, day, lon, lat):
    print(f"Opening ozone file {file}")
    tab_lon = np.arange(-180, 180, 1)
    tab_lat = -np.arange(-90, 90, 1)
    with h5py.File(file) as f:
        data = f["Geophysical Data"][f"ozone_mean_{day:03d}"][()]
    func = RGI((tab_lat, tab_lon), data)
    return func(np.stack([lat, lon], axis=-1))


def get_ancillary(lat, lon, dt: datetime, filename1, filename2, filename3, parmID):
    if filename2 == "":
        return read_climatology(filename1, parmID, dt.strftime("%B"), lat, lon)
    else:
        return read_NRT()


def read_climatology(file1, parm_flag, month, lat, lon):
    tab_lon = np.arange(-180, 180, 1)
    tab_lat = -np.arange(-90, 91, 1)

    def interpolate_data(parameter_key):
        with h5py.File(file1) as f:
            data = f[f"{month}"][parameter_key][()]
            func = RGI((tab_lat, tab_lon), data)
            return func(np.stack([lat, lon], axis=-1))

    match parm_flag:
        case 0:  # wind speed: z
            return interpolate_data("z_wind_mean")
        case 1:  # wind speed: m
            return interpolate_data("m_wind_mean")
        case 2:  # surface pressure
            return interpolate_data("press_mean")
        case 3:  # precipitable water (water vapor)
            return interpolate_data("p_water_mean")
        case 5:  # relative humidity
            return interpolate_data("rel_hum_mean")
        case _:
            print("Unknown parameter ID")
            return None


def read_NRT():
    pass
