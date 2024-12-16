import os

import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d

from aer_struc import AEStr as aestr
from oel_util.libgenutils.genutils import BAD_FLT, BAD_INT
from input_struc import Instr as instr
from l2_struc import L2str
from l12_parms import STDPR

p0 = STDPR


class Geom:
    def __init__(self, solz, senz, phi):
        self.senz = senz
        self.solz = solz
        self.phi = phi
        self.csolz = np.cos(np.radians(self.solz))
        self.csenz = np.cos(np.radians(self.senz))
        self.airmass = 1 / self.csolz + 1 / self.csenz
        self.airmass_plp = None
        self.airmass_sph = None


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

        self.wave = None
        self.solz = None
        self.senz = None
        self.phi = None
        self.scatt = None
        self.dtran_wave = None
        self.dtran_theta = None
        self.dtran_airmass = None

        self.model = []


class Aerosol:
    nsd = 10
    rhtab = []

    @classmethod
    def load_aermod(
        cls,
        sensorID: np.int32,
        wave: np.ndarray,
        nwave: np.int32,
        aermodfile: str,
        models: str,
        nmodels: np.int32,
    ):
        cls.nmodels = nmodels
        print(f"Loading aersol models from {aermodfile}")

        cls.aertab = AerModTabStr(nmodels, sensorID)

        for im in range(nmodels):
            file = os.path.join(f"{aermodfile}_{models[im]}.h5")

            with h5py.File(file) as f:
                cls.aertab.wave = f["wave"][:]
                cls.aertab.solz = f["solz"][:]
                cls.aertab.senz = f["senz"][:]
                cls.aertab.phi = f["phi"][:]
                cls.aertab.scatt = f["scatt"][:]
                cls.aertab.dtran_wave = f["dtran_wave"][:]
                cls.aertab.dtran_theta = f["dtran_theta"][:]

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
                cls.rhtab.append(rh)

            cls.aertab.model.append(aermod)
        cls.rhtab = list(set(cls.rhtab))
        cls.rhtab.sort()
        cls.nrh = len(cls.rhtab)
        print(f"Number of Wavelengths                          {cls.aertab.wave.size}")
        print(f"Number of Solar Zenith Angles                  {cls.aertab.solz.size}")
        print(f"Number of View Zenith Angles                   {cls.aertab.senz.size}")
        print(f"Number of Relative Azimuth Angles              {cls.aertab.phi.size}")
        print(f"Number of Scattering Angles                    {cls.aertab.scatt.size}")
        print(
            f"Number of Diffuse Transmittance Wavelengths    {cls.aertab.dtran_wave.size}"
        )
        print(
            f"Number of Diffuse Transmittance Zenith Angles  {cls.aertab.dtran_theta.size}"
        )
        cls.aertab.dtran_airmass = 1 / np.cos(np.radians(cls.aertab.dtran_theta))

    def ahmader(
        self, sensorID, wave, nwave, iwnir_s, iwnir_l, nmodels, mindx, geom, wv, rhoa
    ):
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
        ) = self.ahmad_atm_corr(
            sensorID, wave, nwave, iwnir_s, iwnir_l, nmodels, mindx, geom, wv, rhoa
        )
        rhoa = rho_aer
        return rhoa, modmin, modmax, modrat, epsnir, tau_pred_min, tau_pred_max

    def rhaer(
        self, sensorID, wave, nwave, iwnir_s, iwnir_l, geom, wv, rh, pr, taur, rhoa
    ):
        # eps1 = 1.0
        # eps2 = 1.0

        # nmodels = self.nsd

        # taua = -np.ones(nwave)
        # tsol = -np.ones(nwave)
        # tsen = -np.ones(nwave)
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
        ) = self.ahmader(
            sensorID, wave, nwave, iwnir_s, iwnir_l, self.nsd, mindx1, geom, wv, rhoa1
        )

        taua1, tsol1, tsen1 = self.diff_tran(
            sensorID,
            wave,
            nwave,
            iwnir_l,
            geom,
            wv,
            pr,
            taur,
            modmin1,
            modmax1,
            modrat1,
            rhoa1,
            tau_pred_min1,
            tau_pred_max1,
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
            ) = self.ahmader(
                sensorID,
                wave,
                nwave,
                iwnir_s,
                iwnir_l,
                self.nsd,
                mindx2,
                geom,
                wv,
                rhoa,
            )

            taua2, tsol2, tsen2 = self.diff_tran(
                sensorID,
                wave,
                nwave,
                iwnir_l,
                geom,
                wv,
                pr,
                taur,
                modmin2,
                modmax2,
                modrat2,
                rhoa2,
                tau_pred_min2,
                tau_pred_max2,
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

    def ahmad_atm_corr(
        self, sensorID, wave, nwave, iwnir_s, iwnir_l, nmodels, mindx, geom, wv, rhoa
    ):

        tau_iwnir_l = np.zeros(nmodels)
        lg_tau_iwnir_s = np.zeros(nmodels)
        tau_iwnir_s = np.zeros(nmodels)
        lg_rho_iwnir_s_pred = np.zeros(nmodels)
        rho_iwnir_s_pred = np.zeros(nmodels)
        eps_pred = np.zeros(nmodels)

        # compute the observed epsilon
        eps_obs = rhoa[iwnir_s] / rhoa[iwnir_l]

        for im in range(nmodels):
            modl = mindx[im]
            # compute AOT at longest aerosol wavelength (iwnir_l)
            ac, bc, cc, dc, ec = self.ms_eps_coef(modl, iwnir_l, wave, geom)

            iwtab_l = iwnir_l
            iwtab_s = iwnir_s

            ax = ac[iwtab_l] - np.log(rhoa[iwnir_l])
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
            eps_pred[im] = rho_iwnir_s_pred[im] / rhoa[iwnir_l]

        epsnir = eps_obs
        im1, im2, mwt = self.model_select_ahmad(nmodels, mindx, eps_pred, eps_obs)
        modmin = mindx[im1]
        modmax = mindx[im2]

        modrat = mwt
        # compute tau_pred and rho_predicted
        tau_pred_min, rho_pred_min = self.comp_rhoa_ms_eps(
            nwave, wave, geom, tau_iwnir_l[im1], modmin, iwnir_l
        )
        tau_pred_max, rho_pred_max = self.comp_rhoa_ms_eps(
            nwave, wave, geom, tau_iwnir_l[im2], modmax, iwnir_l
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

    def comp_rhoa_ms_eps(
        self, nwave: np.int32, wave, geom: Geom, tau_iwnir_l, modl, iwnir_l
    ):
        ac, bc, cc, dc, ec = self.ms_eps_coef(modl, iwnir_l, wave, geom)

        # get the extinction coefficients and compute AOT at all wavelengths

        ext_modl = self.aertab.model[modl].extc
        tau_pred = (ext_modl / ext_modl[iwnir_l]) * tau_iwnir_l
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
        for im in range(self.nmodels):
            if eps_obs < epsilonT[im].eps_obs:
                result_im = im
                break
        im1 = max(min(result_im - 1, self.nmodels - 1), 0)
        im2 = max(min(result_im, self.nmodels - 1), 0)

        # convert table indices to model indices of the input order
        modmin = epsilonT[im1].modnum
        modmax = epsilonT[im2].modnum

        # compute model weighting
        modrat = (eps_obs - epsilonT[im1].eps_obs) / (
            epsilonT[im2].eps_obs - epsilonT[im1].eps_obs
        )

        if modmin == modmax:
            modrat = 1

        return modmin, modmax, modrat

    def ms_eps_coef(self, modnum, iwnir_l, wave, geom: Geom):
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

        solz = np.full_like(wave, geom.solz, np.float64)
        phi = np.full_like(wave, geom.phi, np.float64)
        phi = np.abs(phi)
        senz = np.full_like(wave, geom.senz, np.float64)
        ac = ams_all_func(np.array([wave, solz, phi, senz]).transpose())
        bc = bms_all_func(np.array([wave, solz, phi, senz]).transpose())
        cc = cms_all_func(np.array([wave, solz, phi, senz]).transpose())
        dc = dms_all_func(np.array([wave, solz, phi, senz]).transpose())
        ec = ems_all_func(np.array([wave, solz, phi, senz]).transpose())
        return ac, bc, cc, dc, ec

    def diff_tran(
        self,
        sensorID,
        wave,
        nwave,
        iwnir_l,
        geom,
        wv,
        pr,
        taur,
        modmin,
        modmax,
        modrat,
        rhoa,
        tauamin,
        tauamax,
    ):
        # get diff trans sun to ground and ground to sensor, per band for each model
        tsolmin, tsenmin = self.model_transmittance(
            modmin, wave, nwave, geom.solz, geom.senz, tauamin
        )
        tsolmax, tsenmax = self.model_transmittance(
            modmax, wave, nwave, geom.solz, geom.senz, tauamax
        )
        tsol = tsolmin * (1.0 - modrat) + tsolmax * modrat
        tmp_pressure_diff = np.exp(-0.5 * taur / geom.csolz * (pr / p0 - 1))
        tsol = tsol * tmp_pressure_diff

        tsen = tsenmin * (1.0 - modrat) + tsenmax * modrat
        tmp_pressure_diff = np.exp(-0.5 * taur / geom.csenz * (pr / p0 - 1))
        tsen = tsen * tmp_pressure_diff

        taua = tauamin * (1 - modrat) + tauamax * modrat

        return taua, tsol, tsen

    def model_transmittance(self, modnum, wave, nwave, solz, senz, taua):
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
