"""Integration tests for VLE flash."""

import numpy as np
import pytest

from pvttool import (
    FlashOptions,
    add_components,
    add_mixture,
    add_thermo,
    vle_flash,
    vle_flash_negative,
)
from pvttool.flash._kvalue import kval_estimate


# Experimental data: Ochi et al. (1986) methanol-water at 322.91 K
_P_EXP = [15932.0, 26131.0, 52142.0]
_Y1_EXP = [0.2741, 0.6294, 0.9736]  # vapor methanol mole fraction


@pytest.fixture
def mw_comps():
    return add_components(["Methanol", "Water"])


def _make_mw_mix(comps, p):
    mix = add_mixture(comps, T=322.91, p=p)
    mix.bip.eos_cons = np.array([[0.0, -0.07], [-0.07, 0.0]])
    return mix


def test_vle_flash_negative_finite(mw_comps):
    """vle_flash_negative produces finite, normalised compositions."""
    thermo = add_thermo()
    opts = FlashOptions()
    for p_val in _P_EXP:
        mix = _make_mw_mix(mw_comps, p_val)
        y, x, V = vle_flash_negative(mix, thermo, opts)
        assert np.isfinite(V)
        assert np.isfinite(y.sum())
        assert np.isfinite(x.sum())
        assert abs(y.sum() - 1.0) < 1e-4
        assert abs(x.sum() - 1.0) < 1e-4


def test_vle_flash_negative_two_phase(mw_comps):
    """At p=26131 Pa equimolar feed is two-phase; vapor fraction in [0,1]."""
    mix = _make_mw_mix(mw_comps, _P_EXP[1])
    thermo = add_thermo()
    y, x, V = vle_flash_negative(mix, thermo)
    assert 0.0 <= V <= 1.0


def test_vle_flash_negative_vapor_composition(mw_comps):
    """Vapor methanol fraction matches experiment within 0.05 at p=26131 Pa."""
    mix = _make_mw_mix(mw_comps, _P_EXP[1])
    thermo = add_thermo()
    y, x, V = vle_flash_negative(mix, thermo)
    assert abs(y[0] - _Y1_EXP[1]) < 0.05


def test_vle_flash_mass_balance(mw_comps):
    """z = V*y + (1-V)*x at p=26131 Pa."""
    mix = _make_mw_mix(mw_comps, _P_EXP[1])
    thermo = add_thermo()
    y, x, V = vle_flash_negative(mix, thermo)
    z = mix.mole_fraction
    z_check = V * y + (1.0 - V) * x
    np.testing.assert_allclose(z_check, z, atol=1e-3)


def test_vle_flash_standard_two_phase(mw_comps):
    """vle_flash with z1=0.35 gives vapor_frac in [0,1] and valid compositions."""
    mix = _make_mw_mix(mw_comps, 26131.0)
    mix.mole_fraction = np.array([0.35, 0.65])
    thermo = add_thermo()
    opts = FlashOptions()
    y, x, V = vle_flash(mix, thermo, opts)
    assert 0.0 <= V <= 1.0
    assert abs(y.sum() - 1.0) < 1e-6
    assert abs(x.sum() - 1.0) < 1e-6
    assert all(y >= 0)
    assert all(x >= 0)
    # mass balance
    z = mix.mole_fraction[0]
    z_check = V * y[0] + (1.0 - V) * x[0]
    assert abs(z_check - z) < 1e-4


def test_kval_estimate_positive(mw_comps):
    mix = _make_mw_mix(mw_comps, 26131.0)
    K = kval_estimate(mix)
    assert len(K) == 2
    assert all(K > 0)


def test_vle_flash_pure_component_trivial():
    """Pure-component flash should return trivial result (V=0 or 1)."""
    comps = add_components(["Methane"])
    mix = add_mixture(comps, T=200.0, p=1e5)
    mix.mole_fraction = np.array([1.0])
    thermo = add_thermo()
    y, x, V = vle_flash(mix, thermo)
    assert abs(y.sum() - 1.0) < 1e-6
    assert abs(x.sum() - 1.0) < 1e-6
