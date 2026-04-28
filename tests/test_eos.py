"""Tests for cubic equations of state: PREOS, SRKEOS, PR78EOS."""

import pytest
import numpy as np
from pvt import Component, Mixture, ThermoModel, PREOS, SRKEOS, PR78EOS
from pvt.eos import ResidualProps, select_z_roots


TOL = 1e-10
ATOL = 1e-6


def make_pure(mixture):
    mixture.mole_fraction = np.ones(1)
    return mixture


def make_mixture(names, T, P):
    comps, not_found = Component.from_database_array(names)
    assert len(not_found) == 0, f"Component not found: {not_found}"
    mixture = Mixture(comps, T, P)
    if len(comps) == 1:
        mixture.mole_fraction = np.ones(1)
    return mixture


# ---------------------------------------------------------------------------
# 1. PREOS Z-factors for pure CH4 at subcritical (150 K, 1e6 Pa)
# ---------------------------------------------------------------------------
def test_preos_ch4_subcritical_z():
    mix = make_mixture(['CH4'], 150.0, 1e6)
    thermo = ThermoModel()
    zL, zV, fugacity, HR, props = PREOS(mix, thermo)

    assert zL > 0, f"zL = {zL}"
    assert zV > 0, f"zV = {zV}"
    assert zL < zV, f"zL={zL}, zV={zV}"
    assert np.all(fugacity > 0), f"fugacity={fugacity}"
    assert isinstance(fugacity, np.ndarray)
    assert len(fugacity) == 1


# ---------------------------------------------------------------------------
# 2. PREOS at supercritical (300 K, 10e6 Pa) → single root
# ---------------------------------------------------------------------------
def test_preos_ch4_supercritical_single_root():
    mix = make_mixture(['CH4'], 300.0, 10e6)
    thermo = ThermoModel()
    zL, zV, fugacity, HR, props = PREOS(mix, thermo)

    assert zL > 0
    assert zV > 0
    assert zL == pytest.approx(zV, abs=TOL)


# ---------------------------------------------------------------------------
# 3. Fugacity coefficients should be positive
# ---------------------------------------------------------------------------
def test_fugacity_positive_ch4():
    for T, P in [(150.0, 1e6), (200.0, 5e6), (300.0, 10e6)]:
        mix = make_mixture(['CH4'], T, P)
        thermo = ThermoModel(fugacity_switch=1)
        _zL, _zV, fugacity, _HR, _props = PREOS(mix, thermo)
        assert np.all(fugacity > 0), f"T={T}, P={P} fugacity={fugacity}"


# ---------------------------------------------------------------------------
# 4. SRKEOS subcritical test (150 K, 1e6 Pa)
# ---------------------------------------------------------------------------
def test_srkeos_ch4_subcritical():
    mix = make_mixture(['CH4'], 150.0, 1e6)
    thermo = ThermoModel()
    zL, zV, fugacity, HR, props = SRKEOS(mix, thermo)

    assert zL > 0
    assert zV > 0
    assert zL < zV
    assert np.all(fugacity > 0)


# ---------------------------------------------------------------------------
# 5. PR78EOS should match PREOS for CH4 (omega=0.011 < 0.491)
# ---------------------------------------------------------------------------
def test_pr78eos_matches_preos_ch4():
    mix = make_mixture(['CH4'], 200.0, 3e6)
    thermo = ThermoModel()

    zL_pr, zV_pr, fug_pr, HR_pr, props_pr = PREOS(mix, thermo)
    zL_78, zV_78, fug_78, HR_78, props_78 = PR78EOS(mix, thermo)

    assert zL_pr == pytest.approx(zL_78, rel=1e-12)
    assert zV_pr == pytest.approx(zV_78, rel=1e-12)
    assert HR_pr == pytest.approx(HR_78, rel=1e-12)
    np.testing.assert_allclose(fug_pr, fug_78, rtol=1e-12)


# ---------------------------------------------------------------------------
# 6. Residual enthalpy sign: HR < 0 for liquid
# ---------------------------------------------------------------------------
def test_hr_negative_for_liquid_ch4():
    mix = make_mixture(['CH4'], 150.0, 1e6)
    thermo = ThermoModel(phase=1)  # liquid
    _zL, _zV, _fug, HR, props = PREOS(mix, thermo)
    assert HR < 0, f"Expected HR < 0 for liquid, got HR={HR}"


# ---------------------------------------------------------------------------
# 7. Two-component mixture (CH4/C2H6) works
# ---------------------------------------------------------------------------
def test_two_component_mixture():
    mix = make_mixture(['CH4', 'C2H6'], 250.0, 5e6)
    thermo = ThermoModel()

    zL, zV, fugacity, HR, props = PREOS(mix, thermo)

    assert zL > 0
    assert zV > 0
    assert zL <= zV
    assert len(fugacity) == 2
    assert np.all(np.isfinite(fugacity))
    assert np.all(fugacity > 0)
    assert np.isfinite(HR)


# ---------------------------------------------------------------------------
# 8. compute_props returns ResidualProps namedtuple
# ---------------------------------------------------------------------------
def test_compute_props_returns_residualprops():
    mix = make_mixture(['CH4'], 200.0, 3e6)
    thermo = ThermoModel()
    zL, zV, fugacity, HR, props = PREOS(mix, thermo, compute_props=True)

    assert isinstance(props, ResidualProps)
    assert np.isfinite(props.HR)
    assert np.isfinite(props.SR)
    assert np.isfinite(props.GR)
    assert np.isfinite(props.VR)
    assert np.isfinite(props.Cp_R)
    assert np.isfinite(props.Cv_R)
    assert props.HR == pytest.approx(HR, rel=1e-12)


# ---------------------------------------------------------------------------
# 9. select_z_roots fallback: only complex roots → returns real parts
# ---------------------------------------------------------------------------
def test_select_z_roots_complex_fallback():
    # Three complex-conjugate roots → fallback picks best real part
    roots = np.array([0.1 + 0.5j, 0.1 - 0.5j, 0.8 + 0.01j])
    B = 0.05
    zL, zV = select_z_roots(roots, B)
    assert zL > B
    assert zL == zV, "Single-root fallback when no real physical roots"


# ---------------------------------------------------------------------------
# 10. SRKEOS returns consistent structure
# ---------------------------------------------------------------------------
def test_srkeos_returns_consistent_structure():
    mix = make_mixture(['CH4'], 200.0, 3e6)
    thermo = ThermoModel()
    result = SRKEOS(mix, thermo, compute_props=True)
    assert len(result) == 5
    zL, zV, fugacity, HR, props = result
    assert isinstance(zL, float)
    assert isinstance(zV, float)
    assert isinstance(fugacity, np.ndarray)
    assert isinstance(HR, float)
    assert isinstance(props, ResidualProps)
    assert len(fugacity) == 1
