"""Absolute fugacity helper."""

from __future__ import annotations

import copy

import numpy as np

from pvttool.classes.mixture import Mixture
from pvttool.classes.thermo_model import ThermoModel


def fugacity(
    mixture: Mixture,
    thermo: ThermoModel,
    liquid_x: np.ndarray,
    vapor_y: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute absolute fugacities for liquid and vapor phases.

    Parameters
    ----------
    mixture : Mixture
        Mixture state (temperature, pressure, bip, components).
    thermo : ThermoModel
        EOS configuration.
    liquid_x : np.ndarray, shape (n,)
        Liquid mole fractions.
    vapor_y : np.ndarray, shape (n,)
        Vapor mole fractions.

    Returns
    -------
    liq_fug : np.ndarray, shape (n,)
        Liquid absolute fugacities [Pa].
    vap_fug : np.ndarray, shape (n,)
        Vapor absolute fugacities [Pa].
    """
    p = mixture.pressure
    thermo_liq = copy.copy(thermo)
    thermo_liq.phase = 1
    thermo_liq.fugacity_switch = 1

    thermo_vap = copy.copy(thermo)
    thermo_vap.phase = 2
    thermo_vap.fugacity_switch = 1

    mix_liq = copy.copy(mixture)
    mix_liq.mole_fraction = liquid_x
    _, _, phi_liq, _, _ = thermo_liq.eos(mix_liq, thermo_liq)

    mix_vap = copy.copy(mixture)
    mix_vap.mole_fraction = vapor_y
    _, _, phi_vap, _, _ = thermo_vap.eos(mix_vap, thermo_vap)

    return phi_liq * liquid_x * p, phi_vap * vapor_y * p
