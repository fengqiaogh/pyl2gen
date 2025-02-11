import numpy as np
import h5py
from scipy.interpolate import RegularGridInterpolator as RGI

SOLZNIGHT = 90.0
SOLZNIGHTA = 80.0


class Flag:
    def __init__(self, wave, settings, lon, lat):
        OCDATAROOT = settings["OCDATAROOT"]
        landstr = settings.get("land", "").replace("$OCDATAROOT", OCDATAROOT.as_posix())
        print(f"Loading land mask information from {landstr}")
        with h5py.File(landstr) as f:
            tab_lat = f["lat"][()]
            tab_lon = f["lon"][()]
            watermask = f["watermask"][()]

        func = RGI((tab_lat, tab_lon), watermask)
        water = func((lat, lon)) > 0
        self.land = ~water

        print(f"Loading DEM information from {settings["demfile"]}")
        print(f"Loading ice mask file from {settings["icefile"]}")

        # 计算波长的最小差值，并返回最小差值的索引
        self.ib412 = np.abs(412 - wave).argmin()
        self.ib555 = np.abs(555 - wave).argmin()
        self.ib670 = np.abs(670 - wave).argmin()
        self.ib865 = np.abs(865 - wave).argmin()
        # 计算云层波长的最小差值，并返回最小差值的索引
        self.ibcloud = np.abs(settings["cloud_wave"] - wave).argmin()
        # 打印云层波长的信息
        print(
            f"Using {wave[self.ibcloud]:.1f} nm channel for cloud flagging over water."
        )
        print(f"Using {wave[self.ib412]:.1f} nm channel for cloud flagging over land.")
        # 获取设置中的参数
        self.satzen = settings["satzen"]
        self.sunzen = settings["sunzen"]
        self.glint = settings["glint"]
        self.extreme_glint = settings["extreme_glint"]
        self.albedo = settings["albedo"]

        self.cloud = None
        self.senzmax = None
        self.solzmax = None
        self.mask = None
        self.hilt = None

    def set_mask(self):
        self.mask = np.logical_or(self.land, self.cloud)
        self.mask = np.logical_or(self.mask, self.senzmax)
        self.mask = np.logical_or(self.mask, self.solzmax)
        self.mask = np.logical_or(self.mask, self.hilt)

    def set(self, level1):
        mu0 = level1.csolz

        # Check view angle limits
        self.senzmax = level1.senz > self.satzen
        self.solzmax = level1.solz > self.sunzen

        # 检查是否为非陆地且太阳天顶角小于 SOLZNIGHT
        condition = (~self.land) & (level1.solz < SOLZNIGHT)

        # Check for glint
        glint = condition & (level1.glint_coef > self.glint)
        self.hilt = glint & (level1.glint_coef > self.extreme_glint)

        # Check for clouds (daytime only)
        cloud_albedo = (
            level1.rhos[:, :, self.ibcloud]
            - level1.TLg[:, :, self.ibcloud] * np.pi / mu0 / level1.Fo[self.ibcloud]
        )

        self.cloud = cloud_albedo > self.albedo

        # Set masking
        self.set_mask()
