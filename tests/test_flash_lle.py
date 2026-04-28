"""Tests for liquid-liquid equilibrium flash."""

import numpy as np
import pytest

from pvttool import (
    FlashOptions,
    add_components,
    add_mixture,
    add_thermo,
    lle_flash,
)


def test_lle_flash_returns_valid_compositions():
    """lle_flash returns normalised compositions and fraction in [0,1]."""
    comps = add_components(["Methanol", "Water"])
    mix = add_mixture(comps, T=322.91, p=26131.0)
    mix.mole_fraction = np.array([0.5, 0.5])
    thermo = add_thermo()
    opts = FlashOptions()
    y2, x1, frac = lle_flash(mix, thermo, opts)
    assert 0.0 <= frac <= 1.0
    assert abs(y2.sum() - 1.0) < 1e-5
    assert abs(x1.sum() - 1.0) < 1e-5


def test_lle_flash_mass_balance():
    """Mass balance: z = frac*y2 + (1-frac)*x1."""
    comps = add_components(["Methanol", "Water"])
    mix = add_mixture(comps, T=322.91, p=26131.0)
    mix.mole_fraction = np.array([0.5, 0.5])
    thermo = add_thermo()
    y2, x1, frac = lle_flash(mix, thermo)
    z = mix.mole_fraction
    z_check = frac * y2 + (1.0 - frac) * x1
    np.testing.assert_allclose(z_check, z, atol=1e-3)
