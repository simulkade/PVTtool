"""Shared pytest fixtures."""

import numpy as np
import pytest

from pvttool import (
    Component,
    FlashOptions,
    Mixture,
    ThermoModel,
    add_components,
    add_mixture,
    add_thermo,
    preos,
)


@pytest.fixture
def methane_comp():
    return add_components(["Methane"])[0]


@pytest.fixture
def methane_mix():
    comps = add_components(["Methane"])
    mix = add_mixture(comps, T=300.0, p=1e6)
    mix.mole_fraction = np.array([1.0])
    return mix


@pytest.fixture
def ch4_c10_comps():
    return add_components(["Methane", "n-Decane"])


@pytest.fixture
def ch4_c10_mix(ch4_c10_comps):
    mix = add_mixture(ch4_c10_comps, T=300.0, p=5e6)
    mix.mole_fraction = np.array([0.5, 0.5])
    return mix


@pytest.fixture
def methanol_water_comps():
    return add_components(["Methanol", "Water"])


@pytest.fixture
def methanol_water_mix(methanol_water_comps):
    """Equimolar methanol-water at 322.91 K with kij = -0.07."""
    mix = add_mixture(methanol_water_comps, T=322.91, p=26131.0)
    mix.bip.eos_cons = np.array([[0.0, -0.07], [-0.07, 0.0]])
    return mix


@pytest.fixture
def default_thermo():
    return add_thermo()


@pytest.fixture
def default_options():
    return FlashOptions()
