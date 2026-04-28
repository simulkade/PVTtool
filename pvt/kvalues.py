"""K-value estimation and update functions for flash calculations."""

import numpy as np

from ._constants import R
from .mixture import Mixture
from .thermo_model import ThermoModel


def kval_estimate(mixture: Mixture) -> np.ndarray:
    """Estimate initial K-values using the Wilson correlation.

    K_i = Pc_i / P * exp(5.37 * (1 + omega_i) * (1 - Tc_i / T))

    Args:
        mixture: Mixture object with temperature, pressure, components.

    Returns:
        1D array of estimated K-values.
    """
    Pc = np.array([c.Pc for c in mixture.components])
    Tc = np.array([c.Tc for c in mixture.components])
    omega = np.array([c.acentric_factor for c in mixture.components])
    P = mixture.pressure
    T = mixture.temperature

    K = Pc / P * np.exp(5.37 * (1.0 + omega) * (1.0 - Tc / T))
    return K


def kvalue(
    mixture: Mixture,
    thermo: ThermoModel,
    liquid_x: np.ndarray,
    vapor_y: np.ndarray,
) -> np.ndarray:
    """Compute VLE K-values from EOS fugacity coefficients.

    K_i = phi_L_i / phi_V_i

    Args:
        mixture: Mixture object (used for temperature, pressure, bip, components).
        thermo: ThermoModel object (EOS, mixing rule).
        liquid_x: Liquid mole fractions (1D array).
        vapor_y: Vapor mole fractions (1D array).

    Returns:
        1D array of K-values.
    """
    eos_fn = thermo.EOS

    thermo.phase = 1
    thermo.fugacity_switch = 1
    _, _, liq_fug, _, _ = eos_fn(
        Mixture(mixture.components, mixture.temperature, mixture.pressure,
                liquid_x, mixture.bip),
        thermo, False
    )

    thermo.phase = 2
    _, _, vap_fug, _, _ = eos_fn(
        Mixture(mixture.components, mixture.temperature, mixture.pressure,
                vapor_y, mixture.bip),
        thermo, False
    )

    ki = np.where(vap_fug != 0, liq_fug / vap_fug, 0.0)
    return ki


def kvalueLLE(
    mixture: Mixture,
    thermo: ThermoModel,
    liquid1_x: np.ndarray,
    liquid2_y: np.ndarray,
) -> np.ndarray:
    """Compute LLE K-values from liquid-phase fugacity coefficients.

    K_i = phi_L1_i / phi_L2_i  (both phases evaluated as liquid).

    Args:
        mixture: Mixture object.
        thermo: ThermoModel (phase is set to 1 internally).
        liquid1_x: First liquid phase mole fractions.
        liquid2_y: Second liquid phase mole fractions.

    Returns:
        1D array of LLE K-values.
    """
    eos_fn = thermo.EOS
    thermo.phase = 1
    thermo.fugacity_switch = 1

    _, _, liq1_fug, _, _ = eos_fn(
        Mixture(mixture.components, mixture.temperature, mixture.pressure,
                liquid1_x, mixture.bip),
        thermo, False
    )

    _, _, liq2_fug, _, _ = eos_fn(
        Mixture(mixture.components, mixture.temperature, mixture.pressure,
                liquid2_y, mixture.bip),
        thermo, False
    )

    ki = np.where(liq2_fug != 0, liq1_fug / liq2_fug, 0.0)
    return ki


def fugacity(
    mixture: Mixture,
    thermo: ThermoModel,
    liquid_x: np.ndarray,
    vapor_y: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute partial fugacities (not coefficients) for liquid and vapor phases.

    f_i = phi_i * x_i * P

    Args:
        mixture: Mixture object.
        thermo: ThermoModel.
        liquid_x: Liquid mole fractions.
        vapor_y: Vapor mole fractions.

    Returns:
        Tuple of (liq_fug, vap_fug) — partial fugacity arrays.
    """
    eos_fn = thermo.EOS
    P = mixture.pressure

    thermo.phase = 1
    thermo.fugacity_switch = 1
    _, _, liq_fug_coef, _, _ = eos_fn(
        Mixture(mixture.components, mixture.temperature, mixture.pressure,
                liquid_x, mixture.bip),
        thermo, False
    )

    thermo.phase = 2
    _, _, vap_fug_coef, _, _ = eos_fn(
        Mixture(mixture.components, mixture.temperature, mixture.pressure,
                vapor_y, mixture.bip),
        thermo, False
    )

    return liq_fug_coef * liquid_x * P, vap_fug_coef * vapor_y * P
