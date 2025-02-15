import numpy as np
from libl1.l1 import BANDW


class Bindex:
    WAVE_INDEX_NUM = 13000
    WAVE_INDEX_MIN = 300
    WAVE_INDEX_MAX = WAVE_INDEX_MIN + WAVE_INDEX_NUM

    def __init__(self, wave):
        self.band_index = np.full(self.WAVE_INDEX_NUM + 1, -1, dtype=np.int32)
        self.set(wave)

    def set(self, wave: np.ndarray, dwave_vswir=BANDW):
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
        for wave in [547, 550, 555, 560, 565]:
            ib = self.get(wave)
            if ib >= 0:
                return ib
        return -1


def windex(wave, twave):
    # 设置一个很大的最小差异
    wdiffmin = 99999.0
    index = -1

    # 遍历所有波长
    for iw, tw in enumerate(twave):
        # 如果找到了完全匹配的波长
        if tw == wave:
            index = iw
            break

        # 查找最接近的波长
        wdiff = abs(tw - wave)
        if wdiff < wdiffmin:
            wdiffmin = wdiff
            index = iw

    return index
