"""Peng-Robinson Equation of State (1976)."""

from __future__ import annotations

import numpy as np

from pvttool.classes.mixture import Mixture
from pvttool.classes.thermo_model import ThermoModel
from pvttool.eos._roots import select_z_roots


_R = 8.314  # J/(mol·K)
_S1_PR = 0.623225  # Huron-Vidal constant for PR


def preos(
    mixture: Mixture,
    thermo: ThermoModel,
) -> tuple[float, float, np.ndarray, float, dict]:
    """Peng-Robinson (1976) equation of state.

    Computes compressibility factors, fugacity coefficients, residual
    enthalpy, and optionally a full set of residual properties.

    Parameters
    ----------
    mixture : Mixture
        Mixture state (temperature, pressure, mole_fraction, components, bip).
    thermo : ThermoModel
        Thermodynamic model configuration.

    Returns
    -------
    z_liq : float
        Liquid compressibility factor (smallest real root > B).
    z_vap : float
        Vapor compressibility factor (largest real root > B).
    fugacity : np.ndarray, shape (n,)
        Fugacity coefficients (zeros if ``fugacity_switch == 0``).
    HR : float
        Residual molar enthalpy [J/mol].
    props : dict
        Residual properties with keys ``HR, SR, GR, VR, Cp_R, Cv_R``.
    """
    from pvttool.mixing_rules import mixing_rule

    T = mixture.temperature
    p = mixture.pressure
    x = mixture.mole_fraction
    bip = mixture.bip
    Tc = np.array([c.Tc for c in mixture.components])
    Pc = np.array([c.Pc for c in mixture.components])
    omega = np.array([c.acentric_factor for c in mixture.components])
    N = len(x)

    bi = 0.077796 * _R * Tc / Pc
    aci = 0.457235 * (_R * Tc) ** 2 / Pc
    mi = 0.37646 + (1.54226 - 0.26992 * omega) * omega
    Tr = T / Tc
    alfai = 1.0 + mi * (1.0 - np.sqrt(Tr))
    ai = aci * alfai ** 2

    Q = np.array([-0.53, 0.0]) if thermo.mixing_rule == 3 else (
        np.array([-0.4347, -0.003654]) if thermo.mixing_rule == 4 else np.zeros(2)
    )
    a, b = mixing_rule(mixture, thermo, ai, bi, _S1_PR, Q)

    A_coef = a * p / (_R * T) ** 2
    B_coef = b * p / (_R * T)

    poly = [1.0, -1.0 + B_coef, A_coef - B_coef * (2.0 + 3.0 * B_coef),
            -B_coef * (A_coef - B_coef * (1.0 + B_coef))]
    z_roots = np.roots(poly)
    z_liq, z_vap = select_z_roots(z_roots, B_coef)
    zz = z_liq if thermo.phase == 1 else z_vap

    fugacity = np.zeros(N)
    if thermo.fugacity_switch == 1:
        fugacity = _fugacity_pr(thermo, x, bip, ai, bi, a, b, A_coef, B_coef, zz, T, mixture)

    # da/dT (kij=0 approximation)
    dadT = -(np.sum(x * np.sqrt(ai))) * np.sum(x * np.sqrt(aci) * mi / np.sqrt(Tc)) / np.sqrt(T)
    ln_term = np.log((zz + B_coef * (1.0 + np.sqrt(2.0))) / (zz + B_coef * (1.0 - np.sqrt(2.0))))
    HR = _R * T * (zz - 1.0) + (T * dadT - a) / (b * 2.0 * np.sqrt(2.0)) * ln_term

    props = _residual_props_pr(T, p, x, a, b, ai, aci, mi, Tc, B_coef, zz, dadT, ln_term, HR)

    return z_liq, z_vap, fugacity, HR, props


def _fugacity_pr(thermo, x, bip, ai, bi, a, b, A_coef, B_coef, zz, T, mixture):
    """Compute fugacity coefficients for PR EOS."""
    mr = thermo.mixing_rule
    s1 = _S1_PR
    N = len(x)

    if mr == 1:
        part1 = bi / b * (zz - 1.0) - np.log(zz - B_coef)
        kij = bip.eos_cons + bip.eos_tdep * T
        cross = np.array([
            np.sum(x * np.sqrt(ai[i] * ai) * (1.0 - kij[i, :])) for i in range(N)
        ])
        part3 = (A_coef / (2.828 * B_coef) * (bi / b - 2.0 / a * cross)
                 * np.log((zz + 2.414 * B_coef) / (zz - 0.414 * B_coef)))
        return np.exp(part1 + part3)

    activity_fn = thermo.activity_model
    _, gama = activity_fn(T, x, mixture.components, bip)

    if mr == 2:
        part1 = bi / b * (zz - 1.0) - np.log(zz - B_coef)
        part3 = (-1.0 / (2.0 * np.sqrt(2.0))
                 * (ai / bi / _R / T - np.log(gama) / s1)
                 * np.log((zz + (1.0 + np.sqrt(2.0)) * B_coef) / (zz + (1.0 - np.sqrt(2.0)) * B_coef)))
        return np.exp(part1 + part3)

    alpha = a / (b * _R * T)
    alphai = ai / (bi * _R * T)
    ln_ratio = np.log((zz + (1.0 + np.sqrt(2.0)) * B_coef) / (zz + (1.0 - np.sqrt(2.0)) * B_coef))

    if mr == 3:
        q1 = -0.53
        logfi = (bi / b * (zz - 1.0) - np.log(zz - B_coef)
                 - 1.0 / (2.0 * np.sqrt(2.0))
                 * (alphai + np.log(gama) / q1 + np.log(b / bi) / q1 + (bi / b - 1.0) / q1)
                 * ln_ratio)
        return np.exp(logfi)

    # mr == 4
    q1, q2 = -0.4347, -0.003654
    logfi = (bi / b * (zz - 1.0) - np.log(zz - B_coef)
             - 1.0 / (2.0 * np.sqrt(2.0))
             * (q1 * alphai + q2 * (alpha**2 + alphai**2) + np.log(gama) + np.log(b / bi) + bi / b - 1.0)
             / (q1 + 2.0 * q2 * alpha)
             * ln_ratio)
    return np.exp(logfi)


def _residual_props_pr(T, p, x, a, b, ai, aci, mi, Tc, B_coef, zz, dadT, ln_term, HR):
    """Compute full residual property dict for PR EOS."""
    SR = _R * np.log(zz - B_coef) + dadT / (b * 2.0 * np.sqrt(2.0)) * ln_term
    GR = HR - T * SR
    VR = _R * T * (zz - 1.0) / p

    dsqrtaidT = -np.sqrt(aci) * mi / (2.0 * np.sqrt(Tc * T))
    d2sqrtaidT2 = np.sqrt(aci) * mi / (4.0 * np.sqrt(Tc) * T ** 1.5)
    d2adT2 = 2.0 * np.sum(x * dsqrtaidT) ** 2 + 2.0 * np.sum(x * np.sqrt(ai)) * np.sum(x * d2sqrtaidT2)

    Cv_R = -T * d2adT2 / (b * 2.0 * np.sqrt(2.0)) * ln_term

    V_mol = zz * _R * T / p
    den = V_mol**2 + 2.0 * b * V_mol - b**2
    dPdT_V = _R / (V_mol - b) - dadT / den
    dPdV_T = -_R * T / (V_mol - b) ** 2 + 2.0 * a * (V_mol + b) / den**2
    Cp_R = Cv_R - T * dPdT_V**2 / dPdV_T - _R

    return {"HR": HR, "SR": SR, "GR": GR, "VR": VR, "Cp_R": Cp_R, "Cv_R": Cv_R}
