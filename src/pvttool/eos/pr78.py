"""Peng-Robinson EOS with 1978 alpha-function correction (PR78)."""

from __future__ import annotations

import numpy as np

from pvttool.classes.mixture import Mixture
from pvttool.classes.thermo_model import ThermoModel
from pvttool.eos._roots import select_z_roots
from pvttool.eos.pr import _fugacity_pr, _residual_props_pr


_R = 8.314
_S1_PR = 0.623225


def pr78eos(
    mixture: Mixture,
    thermo: ThermoModel,
) -> tuple[float, float, np.ndarray, float, dict]:
    """Peng-Robinson EOS with 1978 alpha-function correction.

    Identical to PREOS except the m-parameter is modified for components
    with acentric factor > 0.491:

    * omega <= 0.491: m = 0.37646 + 1.54226*omega - 0.26992*omega^2
    * omega > 0.491: m = 0.379642 + 1.48503*omega - 0.164423*omega^2 + 0.016666*omega^3

    Parameters
    ----------
    mixture : Mixture
        Mixture state (temperature, pressure, mole_fraction, components, bip).
    thermo : ThermoModel
        Thermodynamic model configuration.

    Returns
    -------
    z_liq : float
        Liquid compressibility factor.
    z_vap : float
        Vapor compressibility factor.
    fugacity : np.ndarray, shape (n,)
        Fugacity coefficients.
    HR : float
        Residual molar enthalpy [J/mol].
    props : dict
        Residual properties: ``HR, SR, GR, VR, Cp_R, Cv_R``.
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

    bi = 0.0777960739 * _R * Tc / Pc
    aci = 0.457235529 * (_R * Tc) ** 2 / Pc

    # PR78 alpha: different m for omega > 0.491
    mi = np.where(
        omega <= 0.491,
        0.37646 + (1.54226 - 0.26992 * omega) * omega,
        0.379642 + 1.48503 * omega - 0.164423 * omega**2 + 0.016666 * omega**3,
    )

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

    dadT = -(np.sum(x * np.sqrt(ai))) * np.sum(x * np.sqrt(aci) * mi / np.sqrt(Tc)) / np.sqrt(T)
    ln_term = np.log((zz + B_coef * (1.0 + np.sqrt(2.0))) / (zz + B_coef * (1.0 - np.sqrt(2.0))))
    HR = _R * T * (zz - 1.0) + (T * dadT - a) / (b * 2.0 * np.sqrt(2.0)) * ln_term

    props = _residual_props_pr(T, p, x, a, b, ai, aci, mi, Tc, B_coef, zz, dadT, ln_term, HR)

    return z_liq, z_vap, fugacity, HR, props
