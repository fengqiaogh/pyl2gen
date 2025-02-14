import numpy as np


class Bindex:
    def __init__(self):
        self.WAVE_INDEX_NUM = 13000
        self.WAVE_INDEX_MIN = 300
        self.WAVE_INDEX_MAX = self.WAVE_INDEX_MIN + self.WAVE_INDEX_NUM
        self.band_index = np.full(self.WAVE_INDEX_NUM + 1, -1, dtype=np.int32)

    def set(self, wave: np.ndarray, dwave_vswir=10):
        dwave = np.full_like(wave, dwave_vswir, dtype=np.int32)
        dwave[wave > 3000] = 1000
        for iw, wa in enumerate(wave):
            iw1 = (
                np.maximum(wa - dwave[iw] / 2, self.WAVE_INDEX_MIN)
                - self.WAVE_INDEX_MIN
            ).astype(int)
            iw2 = (
                np.minimum(wa + dwave[iw] / 2, self.WAVE_INDEX_MAX)
                - self.WAVE_INDEX_MIN
            ).astype(int)
            if self.band_index[iw1] != -1:
                iw1 = np.maximum(
                    wa - (wa - wave[iw - 1]) / 2,
                    self.WAVE_INDEX_MIN - self.WAVE_INDEX_MIN,
                )
            self.band_index[iw1 : iw2 + 1] = iw

    def get(self, wave: np.int32):
        if wave >= self.WAVE_INDEX_MIN and wave <= self.WAVE_INDEX_MAX:
            return self.band_index[wave - self.WAVE_INDEX_MIN]
        else:
            return -1

    def get_555(self):
        ib = self.get(547)
        if ib < 0:
            ib = self.get(550)
        if ib < 0:
            ib = self.get(555)
        if ib < 0:
            ib = self.get(560)
        if ib < 0:
            ib = self.get(565)

        return ib
