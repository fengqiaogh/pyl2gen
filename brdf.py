import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator as rgi


class Brdf:
    def __init__(self, file):

        with h5py.File(file) as f:
            self.foqtab = f["foq"][()]
            self.phitab = f["phi"][()]
            self.senztab = f["senz"][()]
            self.solztab = f["solz"][()]
            self.chltab = f["chl"][()]
            self.wavetab = f["wave"][()]
            n_a = f["n_phi"][()]
            n_n = f["n_senz"][()]
            n_w = f["n_wave"][()]
            n_c = f["n_chl"][()]
            n_s = f["n_solz"][()]

        print(
            f"morel f/q file dimensions n_a={n_a.size} n_n={n_n.size} n_c={n_c.size} n_s={n_s.size} n_w={n_w.size}"
        )
        print(f"Reading foq file {file}")
        print(f"Closing foq file {file}")
        print(f"Morel f/Q table from file {file}")
        self.lchltab = np.log(self.chltab)

        self.interpolator = rgi(
            (
                self.wavetab,
                self.solztab,
                self.lchltab,
                self.senztab,
                self.phitab,
            ),  # 输入维度
            self.foqtab,  # 查找表数据
            method="linear",  # 线性插值
            bounds_error=False,  # 允许外推
            fill_value=None,  # 外推时使用最近值
        )

    def foqint_morel(self, wave, solz, senzp, phi, chl):
        lchl = np.log(max(chl, 0.01))
        if senzp < self.senztab[0]:
            senzp = self.senztab[0]

        # 构造插值点
        points = np.array(
            [
                wave,
                np.full_like(wave, solz),
                np.full_like(wave, lchl),
                np.full_like(wave, senzp),
                np.full_like(wave, phi),
            ]
        ).T

        brdf = self.interpolator(points)
        return brdf
