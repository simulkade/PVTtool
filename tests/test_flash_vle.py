import numpy as np
import pytest

from pvt import (
    Component,
    FlashOptions,
    Mixture,
    ThermoModel,
    kval_estimate,
    vleflash,
    vleflashnegative,
)


@pytest.fixture
def methanol_water_mixture():
    comps, _ = Component.from_database_array(["CH4O", "H2O"])
    mix = Mixture(comps, 322.91, 26131)
    mix.mole_fraction = np.array([0.5, 0.5])
    mix.bip.EOScons = np.array([[0.0, -0.07], [-0.07, 0.0]])
    return mix


@pytest.fixture
def thermo():
    return ThermoModel()


@pytest.fixture
def opts():
    return FlashOptions()


@pytest.mark.parametrize("pressure", [15932, 26131, 49000])
def test_vleflashnegative_compositions_finite(methanol_water_mixture, thermo, opts, pressure):
    mix = methanol_water_mixture
    mix.pressure = pressure

    y, x, V = vleflashnegative(mix, thermo, opts)

    assert np.all(np.isfinite(y)), f"vapor y not finite at P={pressure}"
    assert np.all(np.isfinite(x)), f"liquid x not finite at P={pressure}"
    assert abs(np.sum(y) - 1.0) < 1e-6, f"vapor y sum != 1 at P={pressure}: {np.sum(y)}"
    assert abs(np.sum(x) - 1.0) < 1e-6, f"liquid x sum != 1 at P={pressure}: {np.sum(x)}"


@pytest.mark.parametrize("pressure", [15932, 26131, 49000])
def test_vleflash_two_phase(methanol_water_mixture, thermo, opts, pressure):
    mix = methanol_water_mixture
    mix.pressure = pressure

    y, x, V = vleflash(mix, thermo, opts)

    assert 0.0 <= V <= 1.0, f"vapor fraction out of range at P={pressure}: {V}"
    assert np.all(np.isfinite(y))
    assert np.all(np.isfinite(x))
    assert abs(np.sum(y) - 1.0) < 1e-6
    assert abs(np.sum(x) - 1.0) < 1e-6


def test_kval_estimate_positive(methanol_water_mixture):
    K = kval_estimate(methanol_water_mixture)

    assert np.all(K > 0), f"K-values not all positive: {K}"
    assert np.all(np.isfinite(K))
