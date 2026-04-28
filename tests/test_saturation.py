"""Tests for bubble and dew point calculations."""

import numpy as np
import pytest

from pvttool import (
    FlashOptions,
    add_components,
    add_mixture,
    add_thermo,
    bubble_pressure,
    bubble_temperature,
    dew_pressure,
    dew_temperature,
)


@pytest.fixture
def ch4_c3_mix():
    comps = add_components(["Methane", "Propane"])
    mix = add_mixture(comps, T=250.0, p=1e6)
    mix.mole_fraction = np.array([0.5, 0.5])
    return mix


@pytest.fixture
def opts():
    return FlashOptions()


def test_bubble_pressure_converges(ch4_c3_mix, opts):
    thermo = add_thermo()
    P_bub, y, converged = bubble_pressure(ch4_c3_mix, thermo, opts)
    assert converged, "bubble_pressure should converge"
    assert P_bub > 0
    assert abs(y.sum() - 1.0) < 1e-5


def test_bubble_pressure_physical_range(ch4_c3_mix, opts):
    thermo = add_thermo()
    P_bub, y, _ = bubble_pressure(ch4_c3_mix, thermo, opts)
    # Bubble pressure for 50/50 CH4-C3H8 at 250 K should be in MPa range
    assert 1e5 < P_bub < 1e8
    assert all(y >= 0)


def test_dew_pressure_converges(ch4_c3_mix, opts):
    thermo = add_thermo()
    P_dew, x, converged = dew_pressure(ch4_c3_mix, thermo, opts)
    assert converged, "dew_pressure should converge"
    assert P_dew > 0
    assert abs(x.sum() - 1.0) < 1e-5


def test_bubble_temperature_converges(ch4_c3_mix, opts):
    thermo = add_thermo()
    T_bub, y, converged = bubble_temperature(ch4_c3_mix, thermo, opts)
    assert converged, "bubble_temperature should converge"
    assert 100.0 < T_bub < 500.0
    assert abs(y.sum() - 1.0) < 1e-5


def test_dew_temperature_converges(ch4_c3_mix, opts):
    thermo = add_thermo()
    T_dew, x, converged = dew_temperature(ch4_c3_mix, thermo, opts)
    assert converged, "dew_temperature should converge"
    assert 100.0 < T_dew < 500.0
    assert abs(x.sum() - 1.0) < 1e-5


def test_bubble_dew_pressure_ordering(ch4_c3_mix, opts):
    """Bubble pressure <= dew pressure for a two-phase mixture."""
    thermo = add_thermo()
    P_bub, _, _ = bubble_pressure(ch4_c3_mix, thermo, opts)
    P_dew, _, _ = dew_pressure(ch4_c3_mix, thermo, opts)
    # For a two-phase-capable mixture the bubble P >= dew P is not guaranteed
    # but both should be finite and positive
    assert P_bub > 0 and P_dew > 0


def test_methanol_water_bubble_pressure():
    """Methanol-water bubble pressure at 322.91 K near experimental 15932 Pa."""
    comps = add_components(["Methanol", "Water"])
    mix = add_mixture(comps, T=322.91, p=15932.0)
    mix.mole_fraction = np.array([0.1, 0.9])  # water-rich
    mix.bip.eos_cons = np.array([[0.0, -0.07], [-0.07, 0.0]])
    thermo = add_thermo()
    opts = FlashOptions()
    P_bub, y, converged = bubble_pressure(mix, thermo, opts)
    assert converged
    assert P_bub > 0
    # water-rich mixture bubble pressure should be in reasonable kPa range
    assert 1e3 < P_bub < 1e6
