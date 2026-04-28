"""Flash calculation and phase stability testing functions.

Provides VLE flash, LLE flash, stability tests, and composition helpers.
All composition vectors are 1D numpy arrays (row vectors from MATLAB perspective).
"""

import numpy as np

from .kvalues import kval_estimate, kvalue, kvalueLLE, fugacity
from .rachford_rice import RachfordRiceNR
from .utils import mynormalize
from .mixture import Mixture
from .thermo_model import ThermoModel
from .flash_options import FlashOptions


def massbalfunc(
    composition: np.ndarray, ki: np.ndarray, vapor_frac: float
) -> tuple[float, float]:
    """Evaluate Rachford-Rice function and its derivative.

    f(V) = sum_i z_i*(K_i-1) / (1 + V*(K_i-1))
    df/dV = -sum_i z_i*(K_i-1)^2 / (1 + V*(K_i-1))^2

    Args:
        composition: [1 x N] mole fractions z_i.
        ki: [1 x N] K-values.
        vapor_frac: Vapor fraction V.

    Returns:
        Tuple of (f, dfdv).
    """
    f = 0.0
    dfdv = 0.0
    for i in range(len(composition)):
        if composition[i] != 0.0:
            denom = 1.0 + vapor_frac * (ki[i] - 1.0)
            f += composition[i] * (ki[i] - 1.0) / denom
            dfdv -= composition[i] * (ki[i] - 1.0) ** 2 / denom ** 2
    return f, dfdv


def xy_calc(
    composition: np.ndarray, vapor_frac: float, ki: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Compute liquid and vapor compositions from feed and K-values.

    liquid_x_i = z_i / (1 + V*(K_i-1))
    vapor_y_i  = K_i * liquid_x_i

    Args:
        composition: Feed mole fractions z_i.
        vapor_frac: Vapor fraction V.
        ki: K-values.

    Returns:
        Tuple of (liquid_x, vapor_y).
    """
    denom = 1.0 + vapor_frac * (ki - 1.0)
    liquid_x = composition / denom
    vapor_y = ki * liquid_x
    return liquid_x, vapor_y


def vleflash(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Isothermal VLE flash at fixed T and P using Rachford-Rice successive substitution.

    Uses GDEM (Dominant Eigenvalue Method) acceleration every 5 iterations.

    Args:
        mixture: Mixture object with T, P, feed composition, BIP.
        thermo: ThermoModel with EOS, mixing rule.
        options: FlashOptions (default if None).

    Returns:
        Tuple of (vapor_y, liquid_x, vapor_frac).
        vapor_y: Vapor mole fractions.
        liquid_x: Liquid mole fractions.
        vapor_frac: Vapor phase fraction [0, 1].
    """
    if options is None:
        options = FlashOptions()

    eps1 = options.accuracy
    max_itr = options.iteration

    ki = kval_estimate(mixture)
    composition = mixture.mole_fraction

    vapor_frac = 0.5
    error1 = 1.0
    error2 = 1.0
    error3 = 1.0

    # GDEM acceleration history
    lnK_prev = None
    lnK_prev2 = None

    j = 0
    while (error1 > eps1 or error2 > eps1 or error3 > eps1):
        j += 1
        if j > max_itr:
            break

        f, dfdv = massbalfunc(composition, ki, vapor_frac)
        vapor_frac -= f / dfdv

        if vapor_frac < 0.0:
            vapor_frac = 0.0
        elif vapor_frac > 1.0:
            vapor_frac = 1.0

        liquid_x, vapor_y = xy_calc(composition, vapor_frac, ki)

        error1 = abs(np.sum(liquid_x) - 1.0)
        error2 = abs(np.sum(vapor_y) - 1.0)
        if abs(dfdv) < eps1:
            error3 = abs(f)
        else:
            error3 = abs(f / dfdv)

        liquid_x = mynormalize(liquid_x)
        vapor_y = mynormalize(vapor_y)

        ki = kvalue(mixture, thermo, liquid_x, vapor_y)

        # GDEM acceleration every 5 iterations
        if j > 0 and j % 5 == 0 and j >= 10:
            lnK_curr = np.log(np.maximum(ki, 1e-20))
            if lnK_prev is not None and lnK_prev2 is not None:
                dg1 = lnK_curr - lnK_prev
                dg2 = lnK_prev - lnK_prev2
                denom = float(np.dot(dg2, dg2))
                if denom > 1e-20:
                    lam = float(np.dot(dg1, dg2)) / denom
                    if 0.01 < lam < 0.99:
                        lnK_acc = lnK_curr + lam / (1.0 - lam) * dg1
                        ki = np.exp(lnK_acc)
            lnK_prev2 = lnK_prev
            lnK_prev = lnK_curr

    if vapor_frac == 0.0 or vapor_frac == 1.0:
        vapor_y = composition.copy()
        liquid_x = composition.copy()

    return mynormalize(vapor_y), mynormalize(liquid_x), vapor_frac


def vleflashnegative(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """VLE flash allowing negative/extended vapor fractions (negative saturation method).

    Used for tracing phase envelopes where feed may be single-phase.
    Converges on fugacity equality rather than mole-fraction sums.

    Args:
        mixture: Mixture object.
        thermo: ThermoModel.
        options: FlashOptions.

    Returns:
        Tuple of (vapor_y, liquid_x, vapor_frac).
    """
    if options is None:
        options = FlashOptions()

    eps1 = options.accuracy
    max_itr = options.iteration

    ki = kval_estimate(mixture)
    composition = mixture.mole_fraction

    # Correct extreme K-values
    Kmax = float(np.max(ki))
    Kmin = float(np.min(ki))
    if Kmax == 0.0:
        Kmax = 1.0
    if Kmax < 1.0:
        idx = np.argmax(ki)
        ki[idx] = 1.0 + Kmax
    if Kmin > 1.0:
        idx = np.argmin(ki)
        ki[idx] = 0.1

    j = 0
    while True:
        j += 1
        if j > max_itr:
            break

        # Physical bounds for vapor fraction
        Kmax = float(np.max(ki))
        Kmin = float(np.min(ki))
        Vmin = -1e10
        Vmax = 1e10
        if Kmax > 1.0:
            Vmin = 1.0 / (1.0 - Kmax)
        if Kmin < 1.0:
            Vmax = 1.0 / (1.0 - Kmin)

        # Solve Rachford-Rice
        vapor_frac, _ = RachfordRiceNR(composition, ki, 0.5)
        if vapor_frac <= Vmin + 1e-10 or vapor_frac >= Vmax - 1e-10:
            # Fallback: bisection
            from .utils import bisection as bisect_fn
            rr_fn = lambda V: massbalfunc(composition, ki, V)[0]
            try:
                vapor_frac = bisect_fn(rr_fn, Vmin + 1e-8, Vmax - 1e-8)
            except RuntimeError:
                pass

        liquid_x, vapor_y = xy_calc(composition, vapor_frac, ki)
        liquid_x = mynormalize(liquid_x)
        vapor_y = mynormalize(vapor_y)

        # Compute fugacities
        liq_fug, vap_fug = fugacity(mixture, thermo, liquid_x, vapor_y)

        # Check convergence on fugacity equality
        ratio = np.where(vap_fug != 0, liq_fug / vap_fug, 1.0)
        error = np.sum((ratio - 1.0) ** 2)
        if error < eps1:
            break

        # Check if single-phase
        if np.sum(np.log(np.maximum(np.abs(ki), 1e-20)) ** 2) < 1e-4:
            if vapor_frac < 0.0 or vapor_frac > 1.0:
                vapor_y = composition.copy()
                liquid_x = composition.copy()
                break

        ki = kvalue(mixture, thermo, liquid_x, vapor_y)

    return vapor_y, liquid_x, vapor_frac


def lleflash(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Liquid-Liquid Equilibrium flash at fixed T and P.

    Both phases are evaluated as liquid. K-values from LLE fugacity ratio.

    Args:
        mixture: Mixture object.
        thermo: ThermoModel.
        options: FlashOptions.

    Returns:
        Tuple of (liquid2_y, liquid1_x, liquid2_frac).
    """
    if options is None:
        options = FlashOptions()

    eps1 = options.accuracy
    max_itr = options.iteration

    ki = kval_estimate(mixture)
    composition = mixture.mole_fraction

    liquid2_frac = 0.5
    error1 = 1.0
    error2 = 1.0
    error3 = 1.0

    j = 0
    while error1 > eps1 or error2 > eps1 or error3 > eps1:
        j += 1
        if j > max_itr:
            break

        f, dfdv = massbalfunc(composition, ki, liquid2_frac)
        liquid2_frac -= f / dfdv

        if liquid2_frac < 0.0:
            liquid2_frac = 0.0
        elif liquid2_frac > 1.0:
            liquid2_frac = 1.0

        liquid1_x, liquid2_y = xy_calc(composition, liquid2_frac, ki)

        error1 = abs(np.sum(liquid1_x) - 1.0)
        error2 = abs(np.sum(liquid2_y) - 1.0)
        if abs(dfdv) < eps1:
            error3 = abs(f)
        else:
            error3 = abs(f / dfdv)

        liquid1_x = mynormalize(liquid1_x)
        liquid2_y = mynormalize(liquid2_y)

        ki = kvalueLLE(mixture, thermo, liquid1_x, liquid2_y)

    if liquid2_frac == 0.0 or liquid2_frac == 1.0:
        liquid2_y = composition.copy()
        liquid1_x = composition.copy()

    return mynormalize(liquid2_y), mynormalize(liquid1_x), liquid2_frac

def stabilityTest(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[np.ndarray, float, float, dict]:
    """Michelsen two-sided stability test for VLE.

    Tests whether a single-phase mixture will spontaneously split into a vapor
    phase (Test 1) or a liquid phase (Test 2).

    Args:
        mixture: Mixture object with T, P, composition, BIP.
        thermo: ThermoModel.
        options: FlashOptions (default if None).

    Returns:
        Tuple of (stability_flag, SL, SV, result).
        stability_flag: 1x2 array; [flag1, flag2]:
            1 = trivial (stable), 2 = non-trivial (unstable), 3 = inconclusive.
        SL: Liquid saturation from Test 2.
        SV: Vapor saturation from Test 1.
        result: Dict with keys 'overall', 'message', 'test1_message', 'test2_message',
                'SL', 'SV'.
    """
    if options is None:
        options = FlashOptions()
    trivial_eps = options.trivialSolutionMaxError
    convergence_eps = options.convergenceMaxError
    max_itr = options.maxIteration

    eosf = thermo.EOS
    thermo.fugacity_switch = 1

    composition = mixture.mole_fraction
    p = mixture.pressure
    T = mixture.temperature

    stability_flag = np.zeros(2, dtype=int)
    SV = 0.0
    SL = 0.0

    # Test 1: can a vapor-like phase form? (feed treated as liquid)
    thermo.phase = 1
    _, _, fug_coef, _, _ = eosf(mixture, thermo)
    mixfug = fug_coef * composition * p

    ki = kval_estimate(mixture)
    conv_error = 1.0 + convergence_eps
    triv_error = 1.0 + trivial_eps
    j = 0
    while conv_error > convergence_eps and triv_error > trivial_eps and j < max_itr:
        j += 1
        Yi = composition * ki
        SV = float(np.sum(Yi))
        yi = Yi / SV
        gas_mix = Mixture(
            mixture.components, mixture.temperature, mixture.pressure, yi, mixture.bip
        )
        thermo.phase = 2
        _, _, fug_coef, _, _ = eosf(gas_mix, thermo)
        gasfug = fug_coef * yi * p
        Ri = mixfug / gasfug * (1.0 / SV)
        ki = ki * Ri
        conv_error = float(np.sum((Ri - 1.0) ** 2))
        triv_error = float(np.sum(np.log(np.maximum(ki, 1e-20)) ** 2))

    if triv_error <= trivial_eps:
        stability_flag[0] = 1
        t1msg = "Test 1 (vapor-like): trivial solution — no vapor phase forms."
    elif conv_error <= convergence_eps:
        stability_flag[0] = 2
        t1msg = "Test 1 (vapor-like): non-trivial convergence — vapor phase can form."
    else:
        stability_flag[0] = 3
        t1msg = "Test 1 (vapor-like): maximum iterations reached — inconclusive."

    # Test 2: can a liquid-like phase form? (feed treated as vapor)
    thermo.phase = 2
    _, _, fug_coef, _, _ = eosf(mixture, thermo)
    mixfug = fug_coef * composition * p

    ki = kval_estimate(mixture)
    conv_error = 1.0 + convergence_eps
    triv_error = 1.0 + trivial_eps
    j = 0
    while conv_error > convergence_eps and triv_error > trivial_eps and j < max_itr:
        j += 1
        Xi = composition / ki
        SL = float(np.sum(Xi))
        xi = Xi / SL
        liq_mix = Mixture(
            mixture.components, mixture.temperature, mixture.pressure, xi, mixture.bip
        )
        thermo.phase = 1
        _, _, fug_coef, _, _ = eosf(liq_mix, thermo)
        liquidfug = fug_coef * xi * p
        Ri = liquidfug / mixfug * SL
        ki = ki * Ri
        conv_error = float(np.sum((Ri - 1.0) ** 2))
        triv_error = float(np.sum(np.log(np.maximum(ki, 1e-20)) ** 2))

    if triv_error <= trivial_eps:
        stability_flag[1] = 1
        t2msg = "Test 2 (liquid-like): trivial solution — no liquid phase forms."
    elif conv_error <= convergence_eps:
        stability_flag[1] = 2
        t2msg = "Test 2 (liquid-like): non-trivial convergence — liquid phase can form."
    else:
        stability_flag[1] = 3
        t2msg = "Test 2 (liquid-like): maximum iterations reached — inconclusive."

    # Build result dict
    if np.any(stability_flag == 2):
        overall = "unstable"
        base = f"Mixture is UNSTABLE. Vapor-liquid split expected at T={T:.2f} K, P={p:.4g} Pa."
    elif np.all(stability_flag == 1):
        overall = "stable"
        base = f"Mixture is stable. No phase split detected at T={T:.2f} K, P={p:.4g} Pa."
    else:
        overall = "inconclusive"
        base = f"Stability result is inconclusive at T={T:.2f} K, P={p:.4g} Pa."

    result = {
        "overall": overall,
        "message": f"{base}\n  {t1msg}\n  {t2msg}",
        "test1_message": t1msg,
        "test2_message": t2msg,
        "SL": SL,
        "SV": SV,
    }

    return stability_flag, SL, SV, result


def stabilityLLETest(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[np.ndarray, float, float, dict]:
    """Michelsen stability test for LLE. Both trial phases evaluated as liquid.

    Args:
        mixture: Mixture object.
        thermo: ThermoModel.
        options: FlashOptions.

    Returns:
        Tuple of (stability_flag, SL, SV, result).
    """
    if options is None:
        options = FlashOptions()
    trivial_eps = options.trivialSolutionMaxError
    convergence_eps = options.convergenceMaxError
    max_itr = options.maxIteration

    eosf = thermo.EOS
    thermo.fugacity_switch = 1
    thermo.phase = 1  # ALL liquid for LLE

    composition = mixture.mole_fraction
    p = mixture.pressure
    T = mixture.temperature

    stability_flag = np.zeros(2, dtype=int)
    SV = 0.0
    SL = 0.0

    def _run_test(y_direction: str) -> tuple[int, float, float, str]:
        """Run one direction of the LLE stability test."""
        nonlocal SV, SL

        thermo.phase = 1
        _, _, fug_coef, _, _ = eosf(mixture, thermo)
        mixfug = fug_coef * composition * p

        ki = kval_estimate(mixture)
        conv_error = 1.0 + convergence_eps
        triv_error = 1.0 + trivial_eps
        j = 0

        while conv_error > convergence_eps and triv_error > trivial_eps and j < max_itr:
            j += 1
            if y_direction == "vapor":
                Yi = composition * ki
                SV = float(np.sum(Yi))
                yi = Yi / SV
                trial_x = yi
                saturation = SV
            else:
                Xi = composition / ki
                SL = float(np.sum(Xi))
                xi = Xi / SL
                trial_x = xi
                saturation = SL

            trial_mix = Mixture(
                mixture.components, mixture.temperature, mixture.pressure,
                trial_x, mixture.bip,
            )
            thermo.phase = 1
            _, _, fug_coef, _, _ = eosf(trial_mix, thermo)
            trialfug = fug_coef * trial_x * p

            if y_direction == "vapor":
                Ri = mixfug / trialfug * (1.0 / saturation)
            else:
                Ri = trialfug / mixfug * saturation
            ki = ki * Ri
            conv_error = float(np.sum((Ri - 1.0) ** 2))
            pos_ki = ki[ki > 0]
            triv_error = float(np.sum(np.log(pos_ki) ** 2)) if len(pos_ki) > 0 else 0.0

        if triv_error <= trivial_eps:
            flag = 1
            msg = f"Test ({y_direction}-like): trivial solution — no {y_direction} phase forms."
        elif conv_error <= convergence_eps:
            flag = 2
            msg = f"Test ({y_direction}-like): non-trivial convergence — {y_direction} phase can form."
        else:
            flag = 3
            msg = f"Test ({y_direction}-like): maximum iterations reached — inconclusive."

        return flag, SV, SL, msg

    stability_flag[0], SV, _, t1msg = _run_test("vapor")
    stability_flag[1], _, SL, t2msg = _run_test("liquid")

    if np.any(stability_flag == 2):
        overall = "unstable"
        base = f"Mixture is UNSTABLE. Liquid-liquid split expected at T={T:.2f} K, P={p:.4g} Pa."
    elif np.all(stability_flag == 1):
        overall = "stable"
        base = f"Mixture is stable. No LLE split detected at T={T:.2f} K, P={p:.4g} Pa."
    else:
        overall = "inconclusive"
        base = f"LLE stability result is inconclusive at T={T:.2f} K, P={p:.4g} Pa."

    result = {
        "overall": overall,
        "message": f"{base}\n  {t1msg}\n  {t2msg}",
        "test1_message": t1msg,
        "test2_message": t2msg,
        "SL": SL,
        "SV": SV,
    }

    return stability_flag, SL, SV, result
