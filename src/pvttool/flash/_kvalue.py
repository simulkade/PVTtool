"""K-value calculations for VLE and LLE."""

from __future__ import annotations

import numpy as np
import copy

from pvttool.classes.mixture import Mixture
from pvttool.classes.thermo_model import ThermoModel


def kval_estimate(mixture: Mixture) -> np.ndarray:
    """Wilson correlation initial K-value estimate.

    Parameters
    ----------
    mixture : Mixture
        Mixture state.

    Returns
    -------
    K : np.ndarray, shape (n,)
        Wilson K-values.
    """
    p = mixture.pressure
    T = mixture.temperature
    Pc = np.array([c.Pc for c in mixture.components])
    Tc = np.array([c.Tc for c in mixture.components])
    omega = np.array([c.acentric_factor for c in mixture.components])
    return (Pc / p) * np.exp(5.37 * (1.0 + omega) * (1.0 - Tc / T))


def kvalue(
    mixture: Mixture,
    thermo: ThermoModel,
    liquid_x: np.ndarray,
    vapor_y: np.ndarray,
) -> np.ndarray:
    """Compute VLE K-values from fugacity coefficients.

    K_i = phi_i^L / phi_i^V where both are computed with ``thermo.eos``.

    Parameters
    ----------
    mixture : Mixture
        Base mixture (temperature, pressure, bip, components).
    thermo : ThermoModel
        EOS configuration.
    liquid_x : np.ndarray, shape (n,)
        Liquid mole fractions.
    vapor_y : np.ndarray, shape (n,)
        Vapor mole fractions.

    Returns
    -------
    K : np.ndarray, shape (n,)
        Equilibrium K-values.
    """
    thermo_liq = copy.copy(thermo)
    thermo_liq.phase = 1
    thermo_liq.fugacity_switch = 1

    thermo_vap = copy.copy(thermo)
    thermo_vap.phase = 2
    thermo_vap.fugacity_switch = 1

    mix_liq = copy.copy(mixture)
    mix_liq.mole_fraction = liquid_x
    _, _, liq_fug, _, _ = thermo_liq.eos(mix_liq, thermo_liq)

    mix_vap = copy.copy(mixture)
    mix_vap.mole_fraction = vapor_y
    _, _, vap_fug, _, _ = thermo_vap.eos(mix_vap, thermo_vap)

    N = len(liquid_x)
    K = np.zeros(N)
    for i in range(N):
        if vap_fug[i] != 0.0:
            K[i] = liq_fug[i] / vap_fug[i]
    return K


def kvalue_lle(
    mixture: Mixture,
    thermo: ThermoModel,
    liquid1_x: np.ndarray,
    liquid2_y: np.ndarray,
) -> np.ndarray:
    """Compute LLE K-values (both phases evaluated as liquid).

    K_i = phi_i^{L1} / phi_i^{L2}

    Parameters
    ----------
    mixture : Mixture
        Base mixture state.
    thermo : ThermoModel
        EOS configuration.
    liquid1_x : np.ndarray, shape (n,)
        Phase 1 mole fractions.
    liquid2_y : np.ndarray, shape (n,)
        Phase 2 mole fractions.

    Returns
    -------
    K : np.ndarray, shape (n,)
        LLE K-values.
    """
    thermo_liq = copy.copy(thermo)
    thermo_liq.phase = 1
    thermo_liq.fugacity_switch = 1

    mix1 = copy.copy(mixture)
    mix1.mole_fraction = liquid1_x
    _, _, liq1_fug, _, _ = thermo_liq.eos(mix1, thermo_liq)

    mix2 = copy.copy(mixture)
    mix2.mole_fraction = liquid2_y
    _, _, liq2_fug, _, _ = thermo_liq.eos(mix2, thermo_liq)

    N = len(liquid1_x)
    K = np.zeros(N)
    for i in range(N):
        if liq2_fug[i] != 0.0:
            K[i] = liq1_fug[i] / liq2_fug[i]
    return K
