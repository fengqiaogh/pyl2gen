import numpy as np


def get_rhos():
    rhos = np.pi / Fo / mu0 * (Lt / tg_sol / tg_sen - Lr) / t_sol / t_sen / t_o2 / t_h2o

    return rhos
