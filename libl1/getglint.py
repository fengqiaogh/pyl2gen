import numpy as np


def getglint_iqu(senz, solz, raz, ws):
    # from Cox & Munk, 54
    ws_gl = 0.049640
    ws_ = np.maximum(ws, 0.001)
    senz_ = np.deg2rad(senz)
    solz_ = np.deg2rad(solz)
    raz_ = np.deg2rad(raz)

    argument = np.cos(senz_) * np.cos(solz_) - np.sin(senz_) * np.sin(solz_) * np.cos(
        raz_
    )
    omega = np.arccos(argument) / 2.0e0
    argument = (np.cos(senz_) + np.cos(solz_)) / (2.0e0 * np.cos(omega))
    beta = np.arccos(argument)
    sigc = ws_gl * np.sqrt(ws_)
    expon = -np.tan(beta) * np.tan(beta) / 2.0 / sigc / sigc
    expon[expon < -30.0] = -30.0  # trap underflow
    expon[expon > 30.0] = 30.0  # trap overflow
    prob = np.exp(expon) / (2.0 * np.pi * sigc * sigc)
    effective_refl, BiRefl = reflec_both(omega)
    cs_beta2 = np.cos(beta) * np.cos(beta)
    cs_beta4 = cs_beta2 * cs_beta2
    glint_coef = effective_refl * prob / (4.0e0 * np.cos(senz_) * cs_beta4)
    CR = (np.cos(solz_) - np.cos(2.0 * omega) * np.cos(senz_)) / (
        np.sin(2.0 * omega) * np.sin(solz_)
    )
    SR = np.sin(solz_) * np.sin(np.pi - raz_) / np.sin(2.0 * omega)
    rot_ang = np.sign(1.0) * np.sign(CR) * np.arcsin(SR)
    c2r = np.cos(2.0 * rot_ang)
    s2r = np.sin(2.0 * rot_ang)
    glint_coef_q = c2r * BiRefl / effective_refl
    glint_coef_u = -s2r * BiRefl / effective_refl
    return glint_coef, glint_coef_q, glint_coef_u


def reflec_both(inc_angle):
    ref = 4.0 / 3.0  # 水的折射率，无波长依赖性
    # 初始化结果矩阵
    Effective_Refl = np.zeros_like(inc_angle)
    BiRefl = np.zeros_like(inc_angle)

    refract_angle = np.arcsin(np.sin(inc_angle) / ref)  # 斯涅尔定律
    Rs = np.sin(inc_angle - refract_angle) / np.sin(inc_angle + refract_angle)
    Rs *= Rs
    Rp = np.tan(inc_angle - refract_angle) / np.tan(inc_angle + refract_angle)
    Rp *= Rp
    Effective_Refl = (Rs + Rp) / 2.0
    BiRefl = (-Rs + Rp) / 2.0

    # 处理 inc_angle < 0.00001 的情况
    small_angle_mask = inc_angle < 0.00001
    Effective_Refl[small_angle_mask] = 0.0204078
    BiRefl[small_angle_mask] = 0.0

    return Effective_Refl, BiRefl
