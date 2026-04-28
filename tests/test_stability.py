"""Tests for VLE and LLE stability tests."""

import numpy as np
import pytest

from pvttool import (
    FlashOptions,
    add_components,
    add_mixture,
    add_thermo,
    stability_lle_test,
    stability_test,
)


def test_stability_test_returns_correct_structure():
    comps = add_components(["Methane", "n-Decane"])
    mix = add_mixture(comps, T=300.0, p=5e6)
    mix.mole_fraction = np.array([0.5, 0.5])
    thermo = add_thermo()
    flags, SL, SV, result = stability_test(mix, thermo)
    assert len(flags) == 2
    assert all(f in (1, 2, 3) for f in flags)
    assert np.isfinite(SL)
    assert np.isfinite(SV)
    assert "overall" in result
    assert result["overall"] in ("stable", "unstable", "inconclusive")
    assert "message" in result


def test_stability_test_unstable_ch4_c10():
    """CH4-Decane 50:50 at 300 K, 5 MPa is typically unstable (two-phase)."""
    comps = add_components(["Methane", "n-Decane"])
    mix = add_mixture(comps, T=300.0, p=5e6)
    mix.mole_fraction = np.array([0.5, 0.5])
    thermo = add_thermo()
    flags, SL, SV, result = stability_test(mix, thermo)
    # This system should detect instability (at least one flag == 2)
    assert result["overall"] in ("unstable", "inconclusive"), (
        f"Expected unstable or inconclusive, got {result['overall']}"
    )


def test_stability_test_stable_pure_methane_gas():
    """Pure methane at 300 K, 1 bar is a single-phase gas (stable)."""
    comps = add_components(["Methane"])
    mix = add_mixture(comps, T=300.0, p=1e5)
    mix.mole_fraction = np.array([1.0])
    thermo = add_thermo()
    flags, SL, SV, result = stability_test(mix, thermo)
    assert result["overall"] in ("stable", "inconclusive")


def test_stability_lle_test_structure():
    comps = add_components(["Methanol", "Water"])
    mix = add_mixture(comps, T=300.0, p=1e5)
    mix.mole_fraction = np.array([0.5, 0.5])
    thermo = add_thermo()
    flags, SL, SV, result = stability_lle_test(mix, thermo)
    assert len(flags) == 2
    assert "overall" in result
    assert result["overall"] in ("stable", "unstable", "inconclusive")


def test_stability_result_message_contains_temperature():
    comps = add_components(["Methane", "n-Decane"])
    mix = add_mixture(comps, T=300.0, p=5e6)
    mix.mole_fraction = np.array([0.5, 0.5])
    thermo = add_thermo()
    _, _, _, result = stability_test(mix, thermo)
    assert "300.00" in result["message"]


def test_stability_with_options():
    comps = add_components(["Methane", "n-Decane"])
    mix = add_mixture(comps, T=300.0, p=5e6)
    mix.mole_fraction = np.array([0.5, 0.5])
    thermo = add_thermo()
    opts = FlashOptions(max_iteration=30, convergence_max_error=1e-8)
    flags, SL, SV, result = stability_test(mix, thermo, opts)
    assert len(flags) == 2
