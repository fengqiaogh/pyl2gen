import numpy as np
from l12_parms import STDPR


def atmocor2(
    wave,
    Lt,
    Lr,
    tLf,
    solz,
    senz,
    phi,
    Tau_r,
    t_sol,
    t_sen,
    tg,
    tg_sol,
    Fo,
    rh,
    pr,
    aerosol,
    fsol,
    Fobar,
    chlorophyll,
    red,
):
    want_ramp = 1
    cbot = 0.7
    ctop = 1.3
    df = 0.33
    p0 = STDPR
    seed_chl = 0.0
    seed_green = 0.0
    seed_red = 0.0
    nir_chg = 0.02

    daer = np.maximum(aerosol.nir_l - aerosol.nir_s, 1)
    cslp = 1.0 / (ctop - cbot)
    cint = -cslp * cbot

    brdf = 1.0

    mu0 = np.cos(np.deg2rad(solz))
    mu = np.cos(np.deg2rad(senz))
    airmass = 1.0 / mu0 + 1.0 / mu

    # Remove pre-computed atmospheric effects

    # Pressure correct the Rayleigh optical thickness
    taur = pr / p0 * Tau_r

    Ltemp = Lt

    # Correct for ozone absorption.  We correct for inbound and outbound here,
    # then we put the inbound back when computing Lw.
    Ltemp = Ltemp / tg

    # Apply polarization correction
    # Ltemp = Ltemp / level1.polcor

    # Remove whitecap radiance
    Ltemp = Ltemp - tLf

    # Subtract the Rayleigh contribution for this geometry.
    Ltemp = Ltemp - Lr

    radref = np.pi / (Fo * mu0)

    # Initialize iteration loop
    chl = seed_chl
    iter = 0
    last_iter = 0
    iter_max = aerosol.aer_iter_max
    iter_min = aerosol.aer_iter_min
    iter_reset = 0
    last_refl_nir = 100.0
    want_glintcorr = 0

    last_tLw_nir = np.zeros_like(Lt)
    rhown_nir = np.zeros_like(Lt)
    tLw_nir = np.zeros_like(Lt)
    Rrs = np.zeros_like(Lt)

    # Begin iterations for aerosol with corrections for non-zero nLw(NIR)
    while not last_iter:
        iter = iter + 1

        # Initialize tLw as surface + aerosol radiance
        tLw = Ltemp

        # Adjust for non-zero NIR water-leaving radiances using IOP model
        rhown_nir = chlorophyll.get_rhown_eval(
            wave,
            Rrs,
            aerosol.nir_s,
            aerosol.nir_l,
            solz,
            senz,
            phi,
            chl,
            rhown_nir,
        )

        # Convert NIR reflectance to TOA W-L radiance
        tLw_nir = rhown_nir / np.pi * Fo * mu0 * t_sol * t_sen / brdf

        # Iteration damping
        tLw_nir = (1.0 - df) * tLw_nir + df * last_tLw_nir

        if chl > 0.0 and chl <= cbot:
            tLw_nir = 0.0
        elif chl > cbot and chl < ctop:
            tLw_nir = tLw_nir * (cslp * chl + cint)

        # Remove estimated NIR water-leaving radiance
        tLw = tLw - tLw_nir

        # Compute the aerosol contribution
        La, t_sol, t_sen, eps, taua = aerosol.get(
            Fo, tLw, mu0, mu, rh, pr, taur, solz, senz, phi
        )

        # Subtract aerosols and normalize
        tLw = tLw - La
        Lw = tLw / t_sen * tg_sol
        nLw = Lw / t_sol / tg_sol / mu0 / fsol * brdf

        # Compute new estimated chlorophyll
        refl_nir = Rrs[red]
        last_tLw_nir[-2:] = tLw_nir[-2:]

        Rrs[:-2] = nLw[:-2] / Fobar[:-2]

        chl = chlorophyll.get(Rrs)

        # Shall we continue iterating
        if iter > iter_max:
            last_iter = 1
        elif (
            np.abs(refl_nir - last_refl_nir) < np.abs(nir_chg * refl_nir)
            or refl_nir < 0.0
        ):
            last_iter = 1
        else:
            pass

        last_refl_nir = refl_nir

    # end of iteration loop
    num_iter = iter

    # Compute final Rrs
    Rrs = nLw / Fobar

    # Compute final chl from final nLw
    chl = chlorophyll.get(Rrs)

    return num_iter, eps, Rrs, chl
