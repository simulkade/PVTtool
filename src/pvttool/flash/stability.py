"""Phase stability tests (Michelsen algorithm)."""

from __future__ import annotations

import copy

import numpy as np

from pvttool.classes.flash_options import FlashOptions
from pvttool.classes.mixture import Mixture
from pvttool.classes.thermo_model import ThermoModel
from pvttool.flash._kvalue import kval_estimate


def stability_test(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[list[int], float, float, dict]:
    """Michelsen two-sided stability test for vapor-liquid equilibrium.

    Tests whether a single-phase mixture will spontaneously split into a
    vapor phase (Test 1) or a liquid phase (Test 2).

    Parameters
    ----------
    mixture : Mixture
        Mixture state (temperature, pressure, mole_fraction, bip, components).
    thermo : ThermoModel
        EOS configuration.
    options : FlashOptions, optional
        Convergence tolerances.  Uses FlashOptions defaults if None.

    Returns
    -------
    stability_flags : list[int]
        Two-element list; each element is:

        * 1 = trivial solution (stable in this direction),
        * 2 = non-trivial convergence (unstable, split expected),
        * 3 = maximum iterations reached (inconclusive).
    SL : float
        Liquid-phase saturation from Test 2.
    SV : float
        Vapor-phase saturation from Test 1.
    result : dict
        Human-readable result with keys ``overall``, ``message``,
        ``test1_message``, ``test2_message``, ``SL``, ``SV``.
    """
    if options is None:
        options = FlashOptions()

    trivial_eps = options.trivial_solution_max_error
    conv_eps = options.convergence_max_error
    max_itr = options.max_iteration

    thermo_liq = copy.copy(thermo)
    thermo_liq.phase = 1
    thermo_liq.fugacity_switch = 1

    thermo_vap = copy.copy(thermo)
    thermo_vap.phase = 2
    thermo_vap.fugacity_switch = 1

    eosf = thermo.eos
    composition = mixture.mole_fraction
    p = mixture.pressure
    stability_flags = [0, 0]
    SV = SL = 0.0

    # -------------------------------------------------------------------------
    # TEST 1: can a vapor-like phase form? (feed treated as liquid)
    # -------------------------------------------------------------------------
    _, _, fug_coef, _, _ = eosf(mixture, thermo_liq)
    mixfug = fug_coef * composition * p

    ki = kval_estimate(mixture)
    conv_error = 1.0 + conv_eps
    triv_error = 1.0 + trivial_eps
    j = 0
    gas_mix = copy.copy(mixture)

    while (conv_error > conv_eps) and (triv_error > trivial_eps) and (j < max_itr):
        j += 1
        Yi = composition * ki
        SV = float(Yi.sum())
        yi = Yi / SV
        gas_mix.mole_fraction = yi
        _, _, fug_coef, _, _ = eosf(gas_mix, thermo_vap)
        gasfug = fug_coef * yi * p
        Ri = mixfug / gasfug * (1.0 / SV)
        ki = ki * Ri
        conv_error = float(np.sum((Ri - 1.0) ** 2))
        pos_ki = ki[ki > 0]
        triv_error = float(np.sum(np.log(pos_ki) ** 2)) if len(pos_ki) else 0.0

    if triv_error <= trivial_eps:
        stability_flags[0] = 1
        t1msg = "Test 1 (vapor-like): trivial solution — no vapor phase forms."
    elif conv_error <= conv_eps:
        stability_flags[0] = 2
        t1msg = "Test 1 (vapor-like): non-trivial convergence — vapor phase can form."
    else:
        stability_flags[0] = 3
        t1msg = "Test 1 (vapor-like): maximum iterations reached — inconclusive."

    # -------------------------------------------------------------------------
    # TEST 2: can a liquid-like phase form? (feed treated as gas)
    # -------------------------------------------------------------------------
    _, _, fug_coef, _, _ = eosf(mixture, thermo_vap)
    mixfug = fug_coef * composition * p

    ki = kval_estimate(mixture)
    conv_error = 1.0 + conv_eps
    triv_error = 1.0 + trivial_eps
    j = 0
    liq_mix = copy.copy(mixture)

    while (conv_error > conv_eps) and (triv_error > trivial_eps) and (j < max_itr):
        j += 1
        Xi = composition / ki
        SL = float(Xi.sum())
        xi = Xi / SL
        liq_mix.mole_fraction = xi
        _, _, fug_coef, _, _ = eosf(liq_mix, thermo_liq)
        liquidfug = fug_coef * xi * p
        Ri = liquidfug / mixfug * SL
        ki = ki * Ri
        conv_error = float(np.sum((Ri - 1.0) ** 2))
        pos_ki = ki[ki > 0]
        triv_error = float(np.sum(np.log(pos_ki) ** 2)) if len(pos_ki) else 0.0

    if triv_error <= trivial_eps:
        stability_flags[1] = 1
        t2msg = "Test 2 (liquid-like): trivial solution — no liquid phase forms."
    elif conv_error <= conv_eps:
        stability_flags[1] = 2
        t2msg = "Test 2 (liquid-like): non-trivial convergence — liquid phase can form."
    else:
        stability_flags[1] = 3
        t2msg = "Test 2 (liquid-like): maximum iterations reached — inconclusive."

    result = _build_result(stability_flags, SL, SV, t1msg, t2msg, mixture)
    return stability_flags, SL, SV, result


def stability_lle_test(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[list[int], float, float, dict]:
    """Michelsen two-sided stability test for liquid-liquid equilibrium.

    Tests whether a single liquid phase will spontaneously split into two
    liquid phases.  Both trial phases are evaluated as liquid (Z-root 1).

    Parameters
    ----------
    mixture : Mixture
        Mixture state.
    thermo : ThermoModel
        EOS configuration.
    options : FlashOptions, optional
        Convergence tolerances.

    Returns
    -------
    stability_flags : list[int]
        Two-element list (1=stable, 2=unstable, 3=inconclusive).
    SL : float
        Phase-B saturation from Test 2.
    SV : float
        Phase-A saturation from Test 1.
    result : dict
        Human-readable result dict.
    """
    if options is None:
        options = FlashOptions()

    trivial_eps = options.trivial_solution_max_error
    conv_eps = options.convergence_max_error
    max_itr = options.max_iteration

    thermo_liq = copy.copy(thermo)
    thermo_liq.phase = 1
    thermo_liq.fugacity_switch = 1

    eosf = thermo.eos
    composition = mixture.mole_fraction
    p = mixture.pressure
    stability_flags = [0, 0]
    SV = SL = 0.0

    # TEST 1
    _, _, fug_coef, _, _ = eosf(mixture, thermo_liq)
    mixfug = fug_coef * composition * p

    ki = kval_estimate(mixture)
    conv_error = 1.0 + conv_eps
    triv_error = 1.0 + trivial_eps
    j = 0
    phase1_mix = copy.copy(mixture)

    while (conv_error > conv_eps) and (triv_error > trivial_eps) and (j < max_itr):
        j += 1
        Yi = composition * ki
        SV = float(Yi.sum())
        yi = Yi / SV
        phase1_mix.mole_fraction = yi
        _, _, fug_coef, _, _ = eosf(phase1_mix, thermo_liq)
        phase1fug = fug_coef * yi * p
        Ri = mixfug / phase1fug * (1.0 / SV)
        ki = ki * Ri
        conv_error = float(np.sum((Ri - 1.0) ** 2))
        pos_ki = ki[ki > 0]
        triv_error = float(np.sum(np.log(pos_ki) ** 2)) if len(pos_ki) else 0.0

    if triv_error <= trivial_eps:
        stability_flags[0] = 1
        t1msg = "Test 1 (phase A): trivial solution — no second liquid phase forms."
    elif conv_error <= conv_eps:
        stability_flags[0] = 2
        t1msg = "Test 1 (phase A): non-trivial convergence — second liquid phase can form."
    else:
        stability_flags[0] = 3
        t1msg = "Test 1 (phase A): maximum iterations reached — inconclusive."

    # TEST 2
    _, _, fug_coef, _, _ = eosf(mixture, thermo_liq)
    mixfug = fug_coef * composition * p

    ki = kval_estimate(mixture)
    conv_error = 1.0 + conv_eps
    triv_error = 1.0 + trivial_eps
    j = 0
    phase2_mix = copy.copy(mixture)

    while (conv_error > conv_eps) and (triv_error > trivial_eps) and (j < max_itr):
        j += 1
        Xi = composition / ki
        SL = float(Xi.sum())
        xi = Xi / SL
        phase2_mix.mole_fraction = xi
        _, _, fug_coef, _, _ = eosf(phase2_mix, thermo_liq)
        phase2fug = fug_coef * xi * p
        Ri = phase2fug / mixfug * SL
        ki = ki * Ri
        conv_error = float(np.sum((Ri - 1.0) ** 2))
        pos_ki = ki[ki > 0]
        triv_error = float(np.sum(np.log(pos_ki) ** 2)) if len(pos_ki) else 0.0

    if triv_error <= trivial_eps:
        stability_flags[1] = 1
        t2msg = "Test 2 (phase B): trivial solution — no second liquid phase forms."
    elif conv_error <= conv_eps:
        stability_flags[1] = 2
        t2msg = "Test 2 (phase B): non-trivial convergence — second liquid phase can form."
    else:
        stability_flags[1] = 3
        t2msg = "Test 2 (phase B): maximum iterations reached — inconclusive."

    result = _build_result(stability_flags, SL, SV, t1msg, t2msg, mixture, lle=True)
    return stability_flags, SL, SV, result


def _build_result(flags, SL, SV, t1msg, t2msg, mixture, lle=False):
    label = "(LLE)" if lle else ""
    T, P = mixture.temperature, mixture.pressure
    if any(f == 2 for f in flags):
        overall = "unstable"
        msg = (f"Mixture is UNSTABLE{' ' + label if label else ''}. "
               f"{'Liquid-liquid' if lle else 'Vapor-liquid'} split expected "
               f"at T={T:.2f} K, P={P:.4g} Pa.\n  {t1msg}\n  {t2msg}")
    elif all(f == 1 for f in flags):
        overall = "stable"
        msg = (f"Mixture is stable{' ' + label if label else ''}. "
               f"No phase split detected at T={T:.2f} K, P={P:.4g} Pa.\n  {t1msg}\n  {t2msg}")
    else:
        overall = "inconclusive"
        msg = (f"Stability result is inconclusive{' ' + label if label else ''} "
               f"at T={T:.2f} K, P={P:.4g} Pa.\n  {t1msg}\n  {t2msg}")

    return {
        "overall": overall,
        "message": msg,
        "test1_message": t1msg,
        "test2_message": t2msg,
        "SL": SL,
        "SV": SV,
    }
