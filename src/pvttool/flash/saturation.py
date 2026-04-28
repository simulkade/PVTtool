"""Bubble and dew point calculations."""

from __future__ import annotations

import copy

import numpy as np

from pvttool.classes.flash_options import FlashOptions
from pvttool.classes.mixture import Mixture
from pvttool.classes.thermo_model import ThermoModel
from pvttool.flash._kvalue import kvalue


def bubble_pressure(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[float, np.ndarray, bool]:
    """Bubble-point pressure at fixed temperature.

    Finds the bubble-point pressure and incipient vapor composition for the
    feed composition at fixed T using successive substitution of EOS K-values.

    Parameters
    ----------
    mixture : Mixture
        Mixture state.  ``mixture.temperature`` is the fixed T [K];
        ``mixture.pressure`` is used as starting neighbourhood.
    thermo : ThermoModel
        EOS configuration.
    options : FlashOptions, optional
        Convergence settings.

    Returns
    -------
    P_bub : float
        Bubble-point pressure [Pa].
    y_bub : np.ndarray, shape (n,)
        Incipient vapor mole fractions.
    converged : bool
        True if the iteration converged.
    """
    if options is None:
        options = FlashOptions()

    T = mixture.temperature
    z = mixture.mole_fraction
    Tc = np.array([c.Tc for c in mixture.components])
    Pc = np.array([c.Pc for c in mixture.components])
    omega = np.array([c.acentric_factor for c in mixture.components])
    eps = options.accuracy

    # Wilson initial P estimate
    P = float(np.sum(z * Pc * np.exp(5.37 * (1.0 + omega) * (1.0 - Tc / T))))
    P = max(P, 1e2)

    K = (Pc / P) * np.exp(5.37 * (1.0 + omega) * (1.0 - Tc / T))
    y = K * z
    S = y.sum()
    y = y / S

    converged = False
    mix1 = copy.copy(mixture)
    for _ in range(options.iteration):
        mix1.pressure = P
        K_new = kvalue(mix1, thermo, z, y)
        S = float(np.dot(K_new, z))
        P = P * S
        P = max(P, 1e2)
        y = K_new * z / S
        y = y / y.sum()

        err_K = float(np.max(np.abs(K_new - K) / np.maximum(np.abs(K), 1e-10)))
        if err_K < eps and abs(S - 1.0) < eps:
            converged = True
            break
        K = K_new

    return P, y, converged


def bubble_temperature(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[float, np.ndarray, bool]:
    """Bubble-point temperature at fixed pressure.

    Finds the bubble-point temperature and incipient vapor composition for
    the feed at fixed P using Wilson Newton pre-conditioning then successive
    substitution of EOS K-values.

    Parameters
    ----------
    mixture : Mixture
        Mixture state.  ``mixture.pressure`` is the fixed P [Pa];
        ``mixture.temperature`` is the starting guess [K].
    thermo : ThermoModel
        EOS configuration.
    options : FlashOptions, optional
        Convergence settings.

    Returns
    -------
    T_bub : float
        Bubble-point temperature [K].
    y_bub : np.ndarray, shape (n,)
        Incipient vapor mole fractions.
    converged : bool
        True if the iteration converged.
    """
    if options is None:
        options = FlashOptions()

    P = mixture.pressure
    z = mixture.mole_fraction
    Tc = np.array([c.Tc for c in mixture.components])
    Pc = np.array([c.Pc for c in mixture.components])
    omega = np.array([c.acentric_factor for c in mixture.components])
    eps = options.accuracy

    # Wilson Newton initial T
    T = mixture.temperature
    for _ in range(30):
        K_w = (Pc / P) * np.exp(5.37 * (1.0 + omega) * (1.0 - Tc / T))
        S = float(np.dot(z, K_w))
        dSdT = float(np.dot(z, K_w * 5.37 * (1.0 + omega) * Tc / T**2))
        dT = (S - 1.0) / max(dSdT, 1e-20)
        T = float(np.clip(T - dT, 50.0, 3000.0))
        if abs(dT) < 0.01:
            break

    K = (Pc / P) * np.exp(5.37 * (1.0 + omega) * (1.0 - Tc / T))
    y = K * z
    y = y / y.sum()

    converged = False
    mix1 = copy.copy(mixture)
    for _ in range(options.iteration):
        mix1.temperature = T
        K_new = kvalue(mix1, thermo, z, y)
        S = float(np.dot(K_new, z))
        y = K_new * z / S
        y = y / y.sum()

        dSdT = float(np.dot(z, K_new * 5.37 * (1.0 + omega) * Tc / T**2))
        T = float(np.clip(T - 0.3 * (S - 1.0) / max(abs(dSdT), 1e-20), 50.0, 3000.0))

        err_K = float(np.max(np.abs(K_new - K) / np.maximum(np.abs(K), 1e-10)))
        if err_K < eps and abs(S - 1.0) < eps:
            converged = True
            break
        K = K_new

    return T, y, converged


def dew_pressure(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[float, np.ndarray, bool]:
    """Dew-point pressure at fixed temperature.

    Parameters
    ----------
    mixture : Mixture
        Mixture state.  ``mixture.temperature`` is the fixed T [K].
    thermo : ThermoModel
        EOS configuration.
    options : FlashOptions, optional
        Convergence settings.

    Returns
    -------
    P_dew : float
        Dew-point pressure [Pa].
    x_dew : np.ndarray, shape (n,)
        Incipient liquid mole fractions.
    converged : bool
        True if the iteration converged.
    """
    if options is None:
        options = FlashOptions()

    T = mixture.temperature
    z = mixture.mole_fraction
    Tc = np.array([c.Tc for c in mixture.components])
    Pc = np.array([c.Pc for c in mixture.components])
    omega = np.array([c.acentric_factor for c in mixture.components])
    eps = options.accuracy

    K_coef = Pc * np.exp(5.37 * (1.0 + omega) * (1.0 - Tc / T))
    P = 1.0 / float(np.sum(z / K_coef))
    P = max(P, 1e2)

    K = K_coef / P
    x = z / K
    x = x / x.sum()

    converged = False
    mix1 = copy.copy(mixture)
    for _ in range(options.iteration):
        mix1.pressure = P
        K_new = kvalue(mix1, thermo, x, z)
        H = float(np.sum(z / K_new))
        P = P / H
        P = max(P, 1e2)
        x = z / K_new
        x = x / x.sum()

        err_K = float(np.max(np.abs(K_new - K) / np.maximum(np.abs(K), 1e-10)))
        if err_K < eps and abs(H - 1.0) < eps:
            converged = True
            break
        K = K_new

    return P, x, converged


def dew_temperature(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[float, np.ndarray, bool]:
    """Dew-point temperature at fixed pressure.

    Parameters
    ----------
    mixture : Mixture
        Mixture state.  ``mixture.pressure`` is the fixed P [Pa];
        ``mixture.temperature`` is the starting guess [K].
    thermo : ThermoModel
        EOS configuration.
    options : FlashOptions, optional
        Convergence settings.

    Returns
    -------
    T_dew : float
        Dew-point temperature [K].
    x_dew : np.ndarray, shape (n,)
        Incipient liquid mole fractions.
    converged : bool
        True if the iteration converged.
    """
    if options is None:
        options = FlashOptions()

    P = mixture.pressure
    z = mixture.mole_fraction
    Tc = np.array([c.Tc for c in mixture.components])
    Pc = np.array([c.Pc for c in mixture.components])
    omega = np.array([c.acentric_factor for c in mixture.components])
    eps = options.accuracy

    # Wilson Newton initial T
    T = mixture.temperature
    for _ in range(30):
        K_w = (Pc / P) * np.exp(5.37 * (1.0 + omega) * (1.0 - Tc / T))
        H = float(np.sum(z / K_w))
        dHdT = -float(np.sum(z / K_w * 5.37 * (1.0 + omega) * Tc / T**2))
        dT = (H - 1.0) / min(dHdT, -1e-20)
        T = float(np.clip(T - dT, 50.0, 3000.0))
        if abs(dT) < 0.01:
            break

    K = (Pc / P) * np.exp(5.37 * (1.0 + omega) * (1.0 - Tc / T))
    x = z / K
    x = x / x.sum()

    converged = False
    mix1 = copy.copy(mixture)
    for _ in range(options.iteration):
        mix1.temperature = T
        K_new = kvalue(mix1, thermo, x, z)
        H = float(np.sum(z / K_new))
        x = z / K_new
        x = x / x.sum()

        dHdT = -float(np.sum(z / K_new * 5.37 * (1.0 + omega) * Tc / T**2))
        T = float(np.clip(T - 0.3 * (H - 1.0) / min(dHdT, -1e-20), 50.0, 3000.0))

        err_K = float(np.max(np.abs(K_new - K) / np.maximum(np.abs(K), 1e-10)))
        if err_K < eps and abs(H - 1.0) < eps:
            converged = True
            break
        K = K_new

    return T, x, converged
