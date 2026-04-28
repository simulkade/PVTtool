"""Vapor-liquid equilibrium flash calculations."""

from __future__ import annotations

import copy

import numpy as np

from pvttool.classes.flash_options import FlashOptions
from pvttool.classes.mixture import Mixture
from pvttool.classes.thermo_model import ThermoModel
from pvttool.flash._kvalue import kval_estimate, kvalue
from pvttool.flash._rachford_rice import mass_bal_func, xy_calc


def _normalize(v: np.ndarray) -> np.ndarray:
    s = v.sum()
    return v / s if s > 0 else v


def vle_flash(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Isothermal vapor-liquid equilibrium flash at fixed T and P.

    Solves the Rachford-Rice equation by Newton-Raphson successive
    substitution with GDEM acceleration every 5 iterations.  Vapor
    fraction is constrained to [0, 1].

    Parameters
    ----------
    mixture : Mixture
        Mixture state (temperature, pressure, mole_fraction, components, bip).
    thermo : ThermoModel
        Thermodynamic model configuration.
    options : FlashOptions, optional
        Convergence tolerances and iteration limits.

    Returns
    -------
    vapor_y : np.ndarray, shape (n,)
        Vapor-phase mole fractions.
    liquid_x : np.ndarray, shape (n,)
        Liquid-phase mole fractions.
    vapor_frac : float
        Molar vapor fraction in [0, 1].
    """
    if options is None:
        options = FlashOptions()

    eps1 = options.accuracy
    max_itr = options.iteration

    ki = kval_estimate(mixture)
    composition = mixture.mole_fraction.copy()

    vapor_frac = 0.5
    error1 = error2 = error3 = 1.0
    j = 0

    lnK_n2 = np.log(np.maximum(ki, 1e-300))
    lnK_n1 = lnK_n2.copy()

    while (error1 > eps1) or (error2 > eps1) or (error3 > eps1):
        j += 1
        if j > max_itr:
            break

        f, dfdv = mass_bal_func(composition, ki, vapor_frac)
        if abs(dfdv) > 0:
            vapor_frac = float(np.clip(vapor_frac - f / dfdv, 0.0, 1.0))

        liquid_x, vapor_y = xy_calc(composition, vapor_frac, ki)

        error1 = abs(liquid_x.sum() - 1.0)
        error2 = abs(vapor_y.sum() - 1.0)
        error3 = f if abs(dfdv) < eps1 else abs(f / dfdv)

        liquid_x = _normalize(liquid_x)
        vapor_y = _normalize(vapor_y)

        ki = kvalue(mixture, thermo, liquid_x, vapor_y)

        # GDEM acceleration every 5 steps (after warm-up)
        lnK_curr = np.log(np.maximum(ki, 1e-300))
        if j % 5 == 0 and j >= 10:
            dg1 = lnK_curr - lnK_n1
            dg2 = lnK_n1 - lnK_n2
            dg2_sq = np.dot(dg2, dg2)
            if dg2_sq > 1e-20:
                lam = np.dot(dg1, dg2) / dg2_sq
                if 0.01 < lam < 0.99:
                    lnK_acc = lnK_curr + lam / (1.0 - lam) * dg1
                    ki = np.exp(lnK_acc)
                    lnK_curr = lnK_acc
        lnK_n2 = lnK_n1.copy()
        lnK_n1 = lnK_curr.copy()

    if vapor_frac == 0.0 or vapor_frac == 1.0:
        vapor_y = composition.copy()
        liquid_x = composition.copy()

    return _normalize(vapor_y), _normalize(liquid_x), vapor_frac


def vle_flash_negative(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """VLE flash allowing vapor fraction outside [0, 1] (negative saturation).

    Identical algorithm to :func:`vle_flash` except the vapor fraction is
    not clamped to [0, 1], allowing extrapolation along the saturation curve
    (useful for bubble/dew point tracing).

    Parameters
    ----------
    mixture : Mixture
        Mixture state.
    thermo : ThermoModel
        Thermodynamic model configuration.
    options : FlashOptions, optional
        Convergence tolerances and iteration limits.

    Returns
    -------
    vapor_y : np.ndarray, shape (n,)
        Vapor-phase mole fractions.
    liquid_x : np.ndarray, shape (n,)
        Liquid-phase mole fractions.
    vapor_frac : float
        Molar vapor fraction (may be outside [0, 1]).
    """
    if options is None:
        options = FlashOptions()

    eps1 = options.accuracy
    max_itr = options.iteration

    ki = kval_estimate(mixture)
    composition = mixture.mole_fraction.copy()

    vapor_frac = 0.5
    error1 = error2 = error3 = 1.0
    j = 0

    lnK_n2 = np.log(np.maximum(ki, 1e-300))
    lnK_n1 = lnK_n2.copy()

    while (error1 > eps1) or (error2 > eps1) or (error3 > eps1):
        j += 1
        if j > max_itr:
            break

        f, dfdv = mass_bal_func(composition, ki, vapor_frac)
        if abs(dfdv) > 0:
            vapor_frac -= f / dfdv

        liquid_x, vapor_y = xy_calc(composition, vapor_frac, ki)

        error1 = abs(liquid_x.sum() - 1.0)
        error2 = abs(vapor_y.sum() - 1.0)
        error3 = f if abs(dfdv) < eps1 else abs(f / dfdv)

        liquid_x = _normalize(liquid_x)
        vapor_y = _normalize(vapor_y)

        ki = kvalue(mixture, thermo, liquid_x, vapor_y)

        lnK_curr = np.log(np.maximum(ki, 1e-300))
        if j % 5 == 0 and j >= 10:
            dg1 = lnK_curr - lnK_n1
            dg2 = lnK_n1 - lnK_n2
            dg2_sq = np.dot(dg2, dg2)
            if dg2_sq > 1e-20:
                lam = np.dot(dg1, dg2) / dg2_sq
                if 0.01 < lam < 0.99:
                    lnK_acc = lnK_curr + lam / (1.0 - lam) * dg1
                    ki = np.exp(lnK_acc)
                    lnK_curr = lnK_acc
        lnK_n2 = lnK_n1.copy()
        lnK_n1 = lnK_curr.copy()

    return _normalize(vapor_y), _normalize(liquid_x), vapor_frac
