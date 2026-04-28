"""Bubble-point and dew-point saturation calculation functions.

Provides bubbleTemperature, bubblePressure, dewTemperature, dewPressure
using successive substitution with EOS K-values.
"""

import numpy as np

from ._constants import R
from .kvalues import kvalue, kval_estimate
from .mixture import Mixture
from .thermo_model import ThermoModel
from .flash_options import FlashOptions
from .utils import mynormalize


def bubbleTemperature(
    mix: Mixture,
    thermo: ThermoModel,
    opts: FlashOptions | None = None,
) -> tuple[float, np.ndarray, int]:
    """Bubble-point temperature at fixed pressure.

    Finds T_bub and incipient vapor composition y_bub for the feed at fixed P.
    Uses Wilson Newton for T initialisation then successive substitution.

    Args:
        mix: Mixture object; mix.pressure is the fixed P, mix.temperature is initial guess.
        thermo: ThermoModel.
        opts: FlashOptions.

    Returns:
        Tuple of (T_bub, y_bub, conv_flag).
        T_bub: Bubble-point temperature [K].
        y_bub: Incipient vapor mole fractions.
        conv_flag: 1 if converged, 0 if not.
    """
    if opts is None:
        opts = FlashOptions()

    composition = mix.mole_fraction
    Tc = np.array([c.Tc for c in mix.components])
    Pc = np.array([c.Pc for c in mix.components])
    omega = np.array([c.acentric_factor for c in mix.components])
    P = mix.pressure
    T = mix.temperature

    # Wilson Newton for initial T
    for _ in range(30):
        Kw = kval_estimate(mix)
        Kw = Pc / P * np.exp(5.37 * (1.0 + omega) * (1.0 - Tc / T))
        S = float(np.dot(composition, Kw))
        dSdT = float(np.dot(composition, Kw * 5.37 * (1.0 + omega) * Tc / T ** 2))
        dT = (S - 1.0) / max(abs(dSdT), 1e-20)
        T = T - dT
        if abs(dT) < 0.01:
            break

    # Successive substitution with T update
    K = Kw
    conv_flag = 0
    for _ in range(opts.iteration):
        S = float(np.dot(composition, K))
        y = K * composition / max(S, 1e-20)

        # Build mixtures for fugacity evaluation
        liq_mix = Mixture(mix.components, T, P, composition, mix.bip)
        vap_mix = Mixture(mix.components, T, P, y, mix.bip)

        K_new = kvalue(liq_mix, thermo, composition, y)

        S_new = float(np.dot(composition, K_new))
        dSdT_new = float(np.dot(composition,
            K_new * 5.37 * (1.0 + omega) * Tc / max(T ** 2, 1e-20)))

        dT = 0.3 * (S_new - 1.0) / max(abs(dSdT_new), 1e-20)
        T = T - dT

        max_err = float(np.max(np.abs(K_new - K) / np.maximum(np.abs(K), 1e-10)))
        K = K_new

        if max_err < opts.accuracy and abs(S_new - 1.0) < opts.accuracy:
            conv_flag = 1
            break

    y_bub = mynormalize(K * composition)
    return float(T), y_bub, conv_flag


def bubblePressure(
    mix: Mixture,
    thermo: ThermoModel,
    opts: FlashOptions | None = None,
) -> tuple[float, np.ndarray, int]:
    """Bubble-point pressure at fixed temperature.

    Finds P_bub and incipient vapor composition y_bub for the feed at fixed T.

    Args:
        mix: Mixture object; mix.temperature is fixed T and mix.pressure is initial guess.
        thermo: ThermoModel.
        opts: FlashOptions.

    Returns:
        Tuple of (P_bub, y_bub, conv_flag).
    """
    if opts is None:
        opts = FlashOptions()

    composition = mix.mole_fraction
    Tc = np.array([c.Tc for c in mix.components])
    Pc = np.array([c.Pc for c in mix.components])
    omega = np.array([c.acentric_factor for c in mix.components])
    T = mix.temperature

    # Wilson initial P
    K_w = Pc / mix.pressure * np.exp(5.37 * (1.0 + omega) * (1.0 - Tc / T))
    P = float(np.dot(composition, K_w * mix.pressure))

    K = Pc / P * np.exp(5.37 * (1.0 + omega) * (1.0 - Tc / T))
    S = float(np.dot(composition, K))
    y = K * composition / max(S, 1e-20)

    conv_flag = 0
    for _ in range(opts.iteration):
        liq_mix = Mixture(mix.components, T, P, composition, mix.bip)
        vap_mix = Mixture(mix.components, T, P, y, mix.bip)

        K_new = kvalue(liq_mix, thermo, composition, y)

        S = float(np.dot(composition, K_new))
        P = P * S
        y = K_new * composition / max(S, 1e-20)

        max_err = float(np.max(np.abs(K_new - K) / np.maximum(np.abs(K), 1e-10)))
        K = K_new

        if max_err < opts.accuracy and abs(S - 1.0) < opts.accuracy:
            conv_flag = 1
            break

    y_bub = mynormalize(K * composition)
    return float(P), y_bub, conv_flag


def dewTemperature(
    mix: Mixture,
    thermo: ThermoModel,
    opts: FlashOptions | None = None,
) -> tuple[float, np.ndarray, int]:
    """Dew-point temperature at fixed pressure.

    Args:
        mix: Mixture object.
        thermo: ThermoModel.
        opts: FlashOptions.

    Returns:
        Tuple of (T_dew, x_dew, conv_flag).
    """
    if opts is None:
        opts = FlashOptions()

    composition = mix.mole_fraction
    Tc = np.array([c.Tc for c in mix.components])
    Pc = np.array([c.Pc for c in mix.components])
    omega = np.array([c.acentric_factor for c in mix.components])
    P = mix.pressure
    T = mix.temperature

    # Wilson Newton for initial T
    for _ in range(30):
        Kw = kval_estimate(mix)
        Kw = Pc / P * np.exp(5.37 * (1.0 + omega) * (1.0 - Tc / T))
        H = float(np.sum(composition / Kw))
        dHdT = -float(np.sum(
            composition / Kw * 5.37 * (1.0 + omega) * Tc / T ** 2
        ))
        dT = (H - 1.0) / max(abs(dHdT), 1e-20)
        T = T - dT
        if abs(dT) < 0.01:
            break

    K = Kw
    conv_flag = 0
    for _ in range(opts.iteration):
        H = float(np.sum(composition / K))
        x = (composition / K) / max(H, 1e-20)

        liq_mix = Mixture(mix.components, T, P, x, mix.bip)
        vap_mix = Mixture(mix.components, T, P, composition, mix.bip)

        K_new = kvalue(liq_mix, thermo, x, composition)

        H_new = float(np.sum(composition / K_new))
        dHdT_new = -float(np.sum(
            composition / K_new * 5.37 * (1.0 + omega) * Tc / max(T, 1e-20) ** 2
        ))

        dT = 0.3 * (H_new - 1.0) / max(abs(dHdT_new), 1e-20)
        T = T - dT

        max_err = float(np.max(np.abs(K_new - K) / np.maximum(np.abs(K), 1e-10)))
        K = K_new

        if max_err < opts.accuracy and abs(H_new - 1.0) < opts.accuracy:
            conv_flag = 1
            break

    x_dew = mynormalize(composition / K)
    return float(T), x_dew, conv_flag


def dewPressure(
    mix: Mixture,
    thermo: ThermoModel,
    opts: FlashOptions | None = None,
) -> tuple[float, np.ndarray, int]:
    """Dew-point pressure at fixed temperature.

    Args:
        mix: Mixture object.
        thermo: ThermoModel.
        opts: FlashOptions.

    Returns:
        Tuple of (P_dew, x_dew, conv_flag).
    """
    if opts is None:
        opts = FlashOptions()

    composition = mix.mole_fraction
    Tc = np.array([c.Tc for c in mix.components])
    Pc = np.array([c.Pc for c in mix.components])
    omega = np.array([c.acentric_factor for c in mix.components])
    T = mix.temperature

    K_w = Pc / mix.pressure * np.exp(5.37 * (1.0 + omega) * (1.0 - Tc / T))
    P = float(np.dot(composition, K_w * mix.pressure))

    K = Pc / P * np.exp(5.37 * (1.0 + omega) * (1.0 - Tc / T))
    H = float(np.sum(composition / K))
    x = (composition / K) / max(H, 1e-20)

    conv_flag = 0
    for _ in range(opts.iteration):
        liq_mix = Mixture(mix.components, T, P, x, mix.bip)
        vap_mix = Mixture(mix.components, T, P, composition, mix.bip)

        K_new = kvalue(liq_mix, thermo, x, composition)

        H = float(np.sum(composition / K_new))
        P = P / H
        x = (composition / K_new) / max(H, 1e-20)

        max_err = float(np.max(np.abs(K_new - K) / np.maximum(np.abs(K), 1e-10)))
        K = K_new

        if max_err < opts.accuracy and abs(H - 1.0) < opts.accuracy:
            conv_flag = 1
            break

    x_dew = mynormalize(composition / K)
    return float(P), x_dew, conv_flag
