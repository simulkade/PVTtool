"""Tests for EOS Z-factors and basic properties."""

import numpy as np
import pytest

from pvttool import (
    add_components,
    add_mixture,
    add_thermo,
    pr78eos,
    preos,
    srkeos,
)
from pvttool.classes.thermo_model import ThermoModel


@pytest.fixture
def pure_methane_mix():
    comps = add_components(["Methane"])
    mix = add_mixture(comps, T=300.0, p=1e6)
    mix.mole_fraction = np.array([1.0])
    return mix


@pytest.fixture
def thermo_no_fug():
    t = add_thermo()
    t.fugacity_switch = 0
    return t


def test_preos_z_factors_pure_methane(pure_methane_mix, thermo_no_fug):
    zl, zv, fug, HR, props = preos(pure_methane_mix, thermo_no_fug)
    assert zv > 0.9, "vapor Z should be near 1 for methane at 300 K, 1 MPa"
    assert zl > 0
    assert zv > zl or zv == zl  # vapor Z >= liquid Z


def test_preos_vapor_z_supercritical_range(pure_methane_mix, thermo_no_fug):
    zl, zv, _, HR, _ = preos(pure_methane_mix, thermo_no_fug)
    # At 300 K, 1 MPa methane is supercritical gas → Z in (0.95, 1.0)
    assert 0.95 < zv < 1.0


def test_srkeos_z_factors(pure_methane_mix, thermo_no_fug):
    t = ThermoModel()
    t.eos = srkeos
    t.fugacity_switch = 0
    zl, zv, fug, HR, props = srkeos(pure_methane_mix, t)
    assert zv > 0.9


def test_pr78eos_matches_preos_low_omega(pure_methane_mix, thermo_no_fug):
    """PR78 should equal PREOS for omega <= 0.491 (methane omega ~ 0.011)."""
    zl_pr, zv_pr, _, _, _ = preos(pure_methane_mix, thermo_no_fug)
    zl_78, zv_78, _, _, _ = pr78eos(pure_methane_mix, thermo_no_fug)
    assert zl_pr == pytest.approx(zl_78, rel=1e-4)
    assert zv_pr == pytest.approx(zv_78, rel=1e-4)


def test_preos_fugacity_coefficients(pure_methane_mix):
    thermo = add_thermo()
    thermo.phase = 2
    _, _, fug, _, _ = preos(pure_methane_mix, thermo)
    assert len(fug) == 1
    assert fug[0] == pytest.approx(0.968, abs=0.02)  # near ideal gas


def test_preos_residual_enthalpy_negative_liquid():
    """Residual enthalpy of liquid should be negative (more ordered than ideal)."""
    comps = add_components(["Propane"])
    mix = add_mixture(comps, T=250.0, p=5e5)  # below critical: liquid
    mix.mole_fraction = np.array([1.0])
    thermo = add_thermo()
    thermo.phase = 1
    thermo.fugacity_switch = 0
    _, _, _, HR, props = preos(mix, thermo)
    assert HR < 0, "Liquid residual enthalpy should be negative"
    assert props["SR"] < 0, "Liquid residual entropy should be negative"
    assert props["VR"] < 0, "Liquid residual volume should be negative"


def test_preos_residual_props_keys():
    comps = add_components(["Methane"])
    mix = add_mixture(comps, T=300.0, p=1e6)
    mix.mole_fraction = np.array([1.0])
    thermo = add_thermo()
    thermo.fugacity_switch = 0
    _, _, _, HR, props = preos(mix, thermo)
    for key in ("HR", "SR", "GR", "VR", "Cp_R", "Cv_R"):
        assert key in props
        assert np.isfinite(props[key])


def test_binary_mixture_preos():
    comps = add_components(["Methane", "Ethane"])
    mix = add_mixture(comps, T=250.0, p=2e6)
    mix.mole_fraction = np.array([0.6, 0.4])
    thermo = add_thermo()
    thermo.phase = 2
    thermo.fugacity_switch = 1
    zl, zv, fug, HR, props = preos(mix, thermo)
    assert len(fug) == 2
    assert all(fug > 0)
    assert all(np.isfinite(fug))
