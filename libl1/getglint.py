import numpy as np


def getglint_iqu(senz, solz, raz, ws):
    # from Cox & Munk, 54
    ws_gl = 0.049640
    senz_ = np.deg2rad(senz)  # 0.640147
