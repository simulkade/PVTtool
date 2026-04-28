"""Tests for residual thermodynamic properties."""

import numpy as np
import pytest

from pvttool import add_components, add_mixture, add_thermo, pr78eos, preos, srkeos
from pvttool.classes.thermo_model import ThermoModel


@pytest.fixture
def liquid_propane_mix():
    comps = add_components(["Propane"])
    mix = add_mixture(comps, T=250.0, p=5e5)
    mix.mole_fraction = np.array([1.0])
    return mix


@pytest.fixture
def thermo_liquid():
    t = add_thermo()
    t.phase = 1
    t.fugacity_switch = 0
    return t


def test_residual_props_all_finite(liquid_propane_mix, thermo_liquid):
    _, _, _, HR, props = preos(liquid_propane_mix, thermo_liquid)
    for key, val in props.items():
        assert np.isfinite(val), f"{key} is not finite: {val}"


def test_hr_negative_for_liquid(liquid_propane_mix, thermo_liquid):
    _, _, _, HR, props = preos(liquid_propane_mix, thermo_liquid)
    assert HR < 0


def test_sr_negative_for_liquid(liquid_propane_mix, thermo_liquid):
    _, _, _, HR, props = preos(liquid_propane_mix, thermo_liquid)
    assert props["SR"] < 0


def test_vr_negative_for_liquid(liquid_propane_mix, thermo_liquid):
    _, _, _, HR, props = preos(liquid_propane_mix, thermo_liquid)
    assert props["VR"] < 0


def test_gr_equals_hr_minus_t_sr(liquid_propane_mix, thermo_liquid):
    T = liquid_propane_mix.temperature
    _, _, _, HR, props = preos(liquid_propane_mix, thermo_liquid)
    assert props["GR"] == pytest.approx(HR - T * props["SR"], rel=1e-6)


def test_hr_prop_equals_fourth_output(liquid_propane_mix, thermo_liquid):
    _, _, _, HR, props = preos(liquid_propane_mix, thermo_liquid)
    assert props["HR"] == pytest.approx(HR, rel=1e-10)


def test_srk_residual_props_finite():
    comps = add_components(["Methane"])
    mix = add_mixture(comps, T=300.0, p=1e6)
    mix.mole_fraction = np.array([1.0])
    t = add_thermo()
    t.eos = srkeos
    t.fugacity_switch = 0
    _, _, _, HR, props = srkeos(mix, t)
    for key, val in props.items():
        assert np.isfinite(val), f"SRK {key} is not finite"


def test_pr78_residual_props_finite():
    comps = add_components(["n-Decane"])
    mix = add_mixture(comps, T=350.0, p=1e5)
    mix.mole_fraction = np.array([1.0])
    t = add_thermo()
    t.eos = pr78eos
    t.phase = 1
    t.fugacity_switch = 0
    _, _, _, HR, props = pr78eos(mix, t)
    assert HR < 0  # liquid decane at 350 K should have negative HR
    for key, val in props.items():
        assert np.isfinite(val)


def test_cp_r_cv_r_sign(liquid_propane_mix, thermo_liquid):
    _, _, _, HR, props = preos(liquid_propane_mix, thermo_liquid)
    # Cp_R and Cv_R for liquid can be positive or negative depending on EOS
    # but must be finite
    assert np.isfinite(props["Cp_R"])
    assert np.isfinite(props["Cv_R"])
