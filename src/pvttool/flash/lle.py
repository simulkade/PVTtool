"""Liquid-liquid equilibrium flash calculation."""

from __future__ import annotations

import numpy as np

from pvttool.classes.flash_options import FlashOptions
from pvttool.classes.mixture import Mixture
from pvttool.classes.thermo_model import ThermoModel
from pvttool.flash._kvalue import kval_estimate, kvalue_lle
from pvttool.flash._rachford_rice import mass_bal_func, xy_calc


def _normalize(v: np.ndarray) -> np.ndarray:
    s = v.sum()
    return v / s if s > 0 else v


def lle_flash(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Liquid-liquid equilibrium flash at fixed T and P.

    Uses successive substitution on the Rachford-Rice equation with LLE
    K-values (K_i = phi_i^{L1} / phi_i^{L2}, both phases evaluated as
    liquid using ``thermo.phase = 1``).

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
    liquid2_y : np.ndarray, shape (n,)
        Phase 2 (second liquid) mole fractions.
    liquid1_x : np.ndarray, shape (n,)
        Phase 1 (first liquid) mole fractions.
    liquid2_frac : float
        Molar fraction of phase 2 in [0, 1].
    """
    if options is None:
        options = FlashOptions()

    eps1 = options.accuracy
    max_itr = options.iteration

    ki = kval_estimate(mixture)
    composition = mixture.mole_fraction.copy()

    liquid2_frac = 0.5
    error1 = error2 = error3 = 1.0
    j = 0

    while (error1 > eps1) or (error2 > eps1) or (error3 > eps1):
        j += 1
        if j > max_itr:
            break

        f, dfdv = mass_bal_func(composition, ki, liquid2_frac)
        liquid2_frac = float(np.clip(liquid2_frac - f / dfdv, 0.0, 1.0))

        liquid1_x, liquid2_y = xy_calc(composition, liquid2_frac, ki)

        error1 = abs(liquid1_x.sum() - 1.0)
        error2 = abs(liquid2_y.sum() - 1.0)
        error3 = f if abs(dfdv) < eps1 else abs(f / dfdv)

        liquid1_x = _normalize(liquid1_x)
        liquid2_y = _normalize(liquid2_y)

        ki = kvalue_lle(mixture, thermo, liquid1_x, liquid2_y)

    if liquid2_frac == 0.0 or liquid2_frac == 1.0:
        liquid2_y = composition.copy()
        liquid1_x = composition.copy()

    return _normalize(liquid2_y), _normalize(liquid1_x), liquid2_frac
