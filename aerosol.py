import os

import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d

from oel_util.libgenutils.genutils import BAD_FLT, BAD_INT
from l12_parms import STDPR

p0 = STDPR


class AerMod:
    def __init__(
        self,
        name,
        rh,
        sd,
        angstrom,
        albedo,
        extc,
        phase,
        acost,
        bcost,
        ccost,
        ams_all,
        bms_all,
        cms_all,
        dms_all,
        ems_all,
        dtran_a,
        dtran_b,
        dtran_a0,
        dtran_b0,
    ):
        self.name = name
        self.rh = rh
        self.sd = sd

        # angstrom exponent (nbands+1)*/
        self.angstrom = angstrom

        # single-scattering albedo(nbands+1), extinction coefficient(nbands+1), phase function/
        self.albedo = albedo
        self.extc = extc
        self.phase = phase

        # quadratic coefficients for SS to MS relationship/
        self.acost = acost
        self.bcost = bcost
        self.ccost = ccost

        # cubic coefficients for ms_epsilon atmospheric correction ..ZA/
        self.ams_all = ams_all
        self.bms_all = bms_all
        self.cms_all = cms_all
        self.dms_all = dms_all
        self.ems_all = ems_all

        # Rayleigh-aerosol diffuse transmittance coeffs/
        self.dtran_a = dtran_a
        self.dtran_b = dtran_b

        self.dtran_a0 = dtran_a0
        self.dtran_b0 = dtran_b0


class AerModTabStr:
    def __init__(self, nmodel, sensorID):
        self.nmodel = nmodel
        self.sensorID = sensorID
        self.model = []


class Aerosol:
    nsd = 10
    rhtab = []
    aer_iter_min = 1

    def __init__(
        self,
        sensorID: np.int32,
        wave: np.ndarray,
        aermodfile: str,
        models: str,
        nmodels: np.int32,
        aer_opt: np.int32,
        iwnir_s_in: np.int32,
        iwnir_l_in: np.int32,
        aer_iter_max: np.int32,
        rhoamin: np.float32,
    ):
        self.sensorID = sensorID
        self.wave = wave
        self.aer_iter_max = aer_iter_max
        self.nmodels = nmodels
        self.nir_s = iwnir_s_in
        self.nir_l = iwnir_l_in
        self.rhoamin = rhoamin
        print(f"Loading aersol models from {aermodfile}")

        self.aertab = AerModTabStr(nmodels, sensorID)

        for im in range(nmodels):
            file = os.path.join(f"{aermodfile}_{models[im]}.h5")

            with h5py.File(file) as f:
                self.aertab.wave = f["wave"][:]
                self.aertab.solz = f["solz"][:]
                self.aertab.senz = f["senz"][:]
                self.aertab.phi = f["phi"][:]
                self.aertab.scatt = f["scatt"][:]
                self.aertab.dtran_wave = f["dtran_wave"][:]
                self.aertab.dtran_theta = f["dtran_theta"][:]

                rh = f.attrs["Relative Humidity"]
                sd = f.attrs["Size Distribution"]
                aermod = AerMod(
                    name=models[im],
                    rh=rh,
                    sd=sd,
                    angstrom=f["angstrom"][:],
                    albedo=f["albedo"][:],
                    extc=f["extc"][:],
                    phase=f["phase"][:],
                    acost=f["acost"][:],
                    bcost=f["bcost"][:],
                    ccost=f["ccost"][:],
                    ams_all=f["ams_all"][:],
                    bms_all=f["bms_all"][:],
                    cms_all=f["cms_all"][:],
                    dms_all=f["dms_all"][:],
                    ems_all=f["ems_all"][:],
                    dtran_a=f["dtran_a"][:],
                    dtran_b=f["dtran_b"][:],
                    dtran_a0=f["dtran_a0"][:],
                    dtran_b0=f["dtran_b0"][:],
                )
            self.rhtab.append(rh)

            self.aertab.model.append(aermod)
        have_rh = 1
        self.rhtab = list(set(self.rhtab))
        self.rhtab.sort()
        self.nrh = len(self.rhtab)
        print(f"Number of Wavelengths                          {self.aertab.wave.size}")
        print(f"Number of Solar Zenith Angles                  {self.aertab.solz.size}")
        print(f"Number of View Zenith Angles                   {self.aertab.senz.size}")
        print(f"Number of Relative Azimuth Angles              {self.aertab.phi.size}")
        print(
            f"Number of Scattering Angles                    {self.aertab.scatt.size}"
        )
        print(
            f"Number of Diffuse Transmittance Wavelengths    {self.aertab.dtran_wave.size}"
        )
        print(
            f"Number of Diffuse Transmittance Zenith Angles  {self.aertab.dtran_theta.size}"
        )
        self.aertab.dtran_airmass = 1 / np.cos(np.radians(self.aertab.dtran_theta))
        if have_rh and self.aertab.nmodel >= 30:
            print("Limiting aerosol models based on RH.")
            use_rh = 1
        match aer_opt:
            case -17:
                print("Using multi-scattering aerosol model selection.")
                print(f"  and NIR correction with up to {self.aer_iter_max} iterations")
                print(
                    f"Using bands at {self.wave[self.nir_s]:4.1f} and {self.wave[self.nir_l]:4.1f} nm for model selection"
                )
                print(f"Extrapolating from {self.wave[self.nir_l]:4.1f} nm")

    def get(self, Fo, La1, csolz, csenz, wv, rh, pr, taur, solz, senz, phi):
        radref = np.pi / Fo / csolz

        # convert input aerosol radiances to relectance
        rhoa = La1 * radref

        # require sufficient signal in two NIR bands
        if (rhoa[-2:] > self.rhoamin).all():
            # require MS epsilon to be reasonable
            eps_tmp = rhoa[self.nir_s] / rhoa[self.nir_l]
            if eps_tmp > 0.1 and eps_tmp < 10:
                (
                    rhoa,
                    modmin1,
                    modmax1,
                    modrat1,
                    modmin2,
                    modmax2,
                    modrat2,
                    eps,
                    taua,
                    t_sol,
                    t_sen,
                ) = self.rhaer(wv, rh, pr, taur, rhoa, solz, senz, phi)
        La2 = rhoa / radref
        return La2, t_sol, t_sen, eps, taua

    def ahmader(self, nmodels, mindx, wv, rhoa, solz, senz, phi):
        (
            modmin,
            modmax,
            modrat,
            epsnir,
            tau_pred_max,
            tau_pred_min,
            rho_pred_max,
            rho_pred_min,
            tau_aer,
            rho_aer,
        ) = self.ahmad_atm_corr(nmodels, mindx, wv, rhoa, solz, senz, phi)
        rhoa = rho_aer
        return rhoa, modmin, modmax, modrat, epsnir, tau_pred_min, tau_pred_max

    def rhaer(self, wv, rh, pr, taur, rhoa, solz, senz, phi):
        rhoa1 = rhoa
        rhoa2 = rhoa
        # rhoa = BAD_FLT

        for i, vrh in enumerate(self.rhtab):
            if vrh > rh:
                irh = i
                break
        irh1 = min(max(0, irh - 1), self.nrh - 2)

        irh2 = irh1 + 1
        wt = (rh - self.rhtab[irh1]) / (self.rhtab[irh2] - self.rhtab[irh1])
        mindx1 = irh1 * self.nsd + np.arange(0, self.nsd)
        mindx2 = irh2 * self.nsd + np.arange(0, self.nsd)

        # compute aerosol reflectances, aot, diffuse trans, eps from first model set
        (
            rhoa1,
            modmin1,
            modmax1,
            modrat1,
            eps1,
            tau_pred_min1,
            tau_pred_max1,
        ) = self.ahmader(self.nsd, mindx1, wv, rhoa1, solz, senz, phi)

        taua1, tsol1, tsen1 = self.diff_tran(
            wv,
            pr,
            taur,
            modmin1,
            modmax1,
            modrat1,
            rhoa1,
            tau_pred_min1,
            tau_pred_max1,
            solz,
            senz,
            phi,
        )

        if irh2 != irh1:

            (
                rhoa2,
                modmin2,
                modmax2,
                modrat2,
                eps2,
                tau_pred_min2,
                tau_pred_max2,
            ) = self.ahmader(self.nsd, mindx2, wv, rhoa, solz, senz, phi)

            taua2, tsol2, tsen2 = self.diff_tran(
                wv,
                pr,
                taur,
                modmin2,
                modmax2,
                modrat2,
                rhoa2,
                tau_pred_min2,
                tau_pred_max2,
                solz,
                senz,
                phi,
            )

            rhoa = (1.0 - wt) * rhoa1 + wt * rhoa2
            taua = (1.0 - wt) * taua1 + wt * taua2
            tsol = (1.0 - wt) * tsol1 + wt * tsol2
            tsen = (1.0 - wt) * tsen1 + wt * tsen2
            eps = (1.0 - wt) * eps1 + wt * eps2
        else:
            rhoa = rhoa1
            taua = taua1
            tsol = tsol1
            tsen = tsen1
            eps = eps1
        return (
            rhoa,
            modmin1,
            modmax1,
            modrat1,
            modmin2,
            modmax2,
            modrat2,
            eps,
            taua,
            tsol,
            tsen,
        )

    def ahmad_atm_corr(self, nmodels, mindx, wv, rhoa, solz, senz, phi):

        tau_iwnir_l = np.zeros(nmodels)
        lg_tau_iwnir_s = np.zeros(nmodels)
        tau_iwnir_s = np.zeros(nmodels)
        lg_rho_iwnir_s_pred = np.zeros(nmodels)
        rho_iwnir_s_pred = np.zeros(nmodels)
        eps_pred = np.zeros(nmodels)

        # compute the observed epsilon
        eps_obs = rhoa[self.nir_s] / rhoa[self.nir_l]

        for im in range(nmodels):
            modl = mindx[im]
            # compute AOT at longest aerosol wavelength (iwnir_l)
            ac, bc, cc, dc, ec = self.ms_eps_coef(modl, solz, senz, phi)

            iwtab_l = self.nir_l
            iwtab_s = self.nir_s

            ax = ac[iwtab_l] - np.log(rhoa[self.nir_l])
            bx = bc[iwtab_l]
            cx = cc[iwtab_l]
            fx = bx * bx - 4.0 * ax * cx
            lg_tau_iwnir_l = 0.5 * (-bx + np.sqrt(fx)) / cx
            tau_iwnir_l[im] = np.exp(lg_tau_iwnir_l)

            # compute AOT at shortest aerosol wavelength (iwnir_s)
            ext_iwnir_l = self.aertab.model[modl].extc[iwtab_l]
            ext_iwnir_s = self.aertab.model[modl].extc[iwtab_s]
            tau_iwnir_s[im] = (ext_iwnir_s / ext_iwnir_l) * tau_iwnir_l[im]
            lg_tau_iwnir_s[im] = np.log(tau_iwnir_s[im])

            # compute reflectance at (iwnir_s)
            lg_rho_iwnir_s_pred[im] = (
                ac[iwtab_s]
                + bc[iwtab_s] * lg_tau_iwnir_s[im]
                + cc[iwtab_s] * np.power(lg_tau_iwnir_s[im], 2)
                + dc[iwtab_s] * np.power(lg_tau_iwnir_s[im], 3)
                + ec[iwtab_s] * np.power(lg_tau_iwnir_s[im], 4)
            )
            rho_iwnir_s_pred[im] = np.exp(lg_rho_iwnir_s_pred[im])

            # compute model epsilon
            eps_pred[im] = rho_iwnir_s_pred[im] / rhoa[self.nir_l]

        epsnir = eps_obs
        im1, im2, mwt = self.model_select_ahmad(nmodels, mindx, eps_pred, eps_obs)
        modmin = mindx[im1]
        modmax = mindx[im2]

        modrat = mwt
        # compute tau_pred and rho_predicted
        tau_pred_min, rho_pred_min = self.comp_rhoa_ms_eps(
            tau_iwnir_l[im1], modmin, solz, senz, phi
        )
        tau_pred_max, rho_pred_max = self.comp_rhoa_ms_eps(
            tau_iwnir_l[im2], modmax, solz, senz, phi
        )

        # compute weighted tau_aer and rho_aer
        tau_aer = (1.0 - mwt) * tau_pred_min + mwt * tau_pred_max
        rho_aer = (1.0 - mwt) * rho_pred_min + mwt * rho_pred_max

        return (
            modmin,
            modmax,
            modrat,
            epsnir,
            tau_pred_max,
            tau_pred_min,
            rho_pred_max,
            rho_pred_min,
            tau_aer,
            rho_aer,
        )

    def comp_rhoa_ms_eps(self, tau_iwnir_l, modl, solz, senz, phi):
        ac, bc, cc, dc, ec = self.ms_eps_coef(modl, solz, senz, phi)

        # get the extinction coefficients and compute AOT at all wavelengths

        ext_modl = self.aertab.model[modl].extc
        tau_pred = (ext_modl / ext_modl[self.nir_l]) * tau_iwnir_l
        lg_tau_pred = np.log(tau_pred)

        # compute rho_pred
        lg_rho_pred = (
            ac
            + bc * lg_tau_pred
            + cc * np.power(lg_tau_pred, 2)
            + dc * np.power(lg_tau_pred, 3)
            + ec * np.power(lg_tau_pred, 4)
        )

        rho_pred = np.exp(lg_rho_pred)
        return tau_pred, rho_pred

    def model_select_ahmad(self, nmodels: np.int32, mindx: np.int32, eps_pred, eps_obs):
        epsilonT = [EpsilonT(modnum=im, eps_obs=eps_pred[im]) for im in range(nmodels)]

        # sort into ascending model epsilon order
        epsilonT = sorted(epsilonT, key=lambda x: x.eps_obs)

        # find bounding epsilon indices in table
        for im in range(nmodels):
            if eps_obs < epsilonT[im].eps_obs:
                break
            else:
                im = nmodels
        im1 = max(min(im - 1, nmodels - 1), 0)
        im2 = max(min(im, nmodels - 1), 0)

        # convert table indices to model indices of the input order
        modmin = epsilonT[im1].modnum
        modmax = epsilonT[im2].modnum

        # compute model weighting
        if modmin == modmax:
            modrat = 1
        else:
            modrat = (eps_obs - epsilonT[im1].eps_obs) / (
                epsilonT[im2].eps_obs - epsilonT[im1].eps_obs
            )

        return modmin, modmax, modrat

    def ms_eps_coef(self, modnum, solz_in, senz_in, phi_in):
        lut_wave = self.aertab.wave
        lut_solz = self.aertab.solz
        lut_phi = self.aertab.phi
        lut_senz = self.aertab.senz

        ams_all = self.aertab.model[modnum].ams_all
        bms_all = self.aertab.model[modnum].bms_all
        cms_all = self.aertab.model[modnum].cms_all
        dms_all = self.aertab.model[modnum].dms_all
        ems_all = self.aertab.model[modnum].ems_all

        ams_all_func = RegularGridInterpolator(
            (lut_wave, lut_solz, lut_phi, lut_senz), ams_all
        )
        bms_all_func = RegularGridInterpolator(
            (lut_wave, lut_solz, lut_phi, lut_senz), bms_all
        )
        cms_all_func = RegularGridInterpolator(
            (lut_wave, lut_solz, lut_phi, lut_senz), cms_all
        )
        dms_all_func = RegularGridInterpolator(
            (lut_wave, lut_solz, lut_phi, lut_senz), dms_all
        )
        ems_all_func = RegularGridInterpolator(
            (lut_wave, lut_solz, lut_phi, lut_senz), ems_all
        )

        solz = np.full_like(self.wave, solz_in, np.float64)
        phi = np.full_like(self.wave, phi_in, np.float64)
        phi = np.abs(phi)
        senz = np.full_like(self.wave, senz_in, np.float64)
        ac = ams_all_func(np.array([self.wave, solz, phi, senz]).transpose())
        bc = bms_all_func(np.array([self.wave, solz, phi, senz]).transpose())
        cc = cms_all_func(np.array([self.wave, solz, phi, senz]).transpose())
        dc = dms_all_func(np.array([self.wave, solz, phi, senz]).transpose())
        ec = ems_all_func(np.array([self.wave, solz, phi, senz]).transpose())
        return ac, bc, cc, dc, ec

    def diff_tran(
        self,
        wv,
        pr,
        taur,
        modmin,
        modmax,
        modrat,
        rhoa,
        tauamin,
        tauamax,
        solz,
        senz,
        phi,
    ):
        csolz = np.cos(np.deg2rad(solz))
        csenz = np.cos(np.deg2rad(senz))

        # get diff trans sun to ground and ground to sensor, per band for each model
        tsolmin, tsenmin = self.model_transmittance(modmin, solz, senz, tauamin)
        tsolmax, tsenmax = self.model_transmittance(modmax, solz, senz, tauamax)
        tsol = tsolmin * (1.0 - modrat) + tsolmax * modrat
        tmp_pressure_diff = np.exp(-0.5 * taur / csolz * (pr / p0 - 1))
        tsol = tsol * tmp_pressure_diff

        tsen = tsenmin * (1.0 - modrat) + tsenmax * modrat
        tmp_pressure_diff = np.exp(-0.5 * taur / csenz * (pr / p0 - 1))
        tsen = tsen * tmp_pressure_diff

        taua = tauamin * (1 - modrat) + tauamax * modrat

        return taua, tsol, tsen

    def model_transmittance(self, modnum, solz, senz, taua):
        dtran_theta = self.aertab.dtran_theta
        dtran_a0 = self.aertab.model[modnum].dtran_a0
        dtran_b0 = self.aertab.model[modnum].dtran_b0
        dtran_a = self.aertab.model[modnum].dtran_a
        dtran_b = self.aertab.model[modnum].dtran_b

        a0_f = interp1d(dtran_theta, dtran_a0)
        b0_f = interp1d(dtran_theta, dtran_b0)
        a_f = interp1d(dtran_theta, dtran_a)
        b_f = interp1d(dtran_theta, dtran_b)

        dtran0 = a0_f(solz) * np.exp(-b0_f(solz) * taua)
        dtran = a_f(senz) * np.exp(-b_f(senz) * taua)
        return dtran0, dtran


class EpsilonT:
    def __init__(self, modnum, eps_obs):
        self.modnum = modnum
        self.eps_obs = eps_obs
