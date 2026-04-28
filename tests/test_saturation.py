import numpy as np
import pytest

from pvt import (
    Component,
    FlashOptions,
    Mixture,
    ThermoModel,
    bubblePressure,
    bubbleTemperature,
    dewPressure,
    dewTemperature,
)


@pytest.fixture
def thermo():
    return ThermoModel()


@pytest.fixture
def opts():
    return FlashOptions()


@pytest.fixture
def methanol_water_mixture():
    comps, _ = Component.from_database_array(["CH4O", "H2O"])
    mix = Mixture(comps, 350.0, 1e5)
    mix.mole_fraction = np.array([0.5, 0.5])
    mix.bip.EOScons = np.array([[0.0, -0.07], [-0.07, 0.0]])
    return mix


def test_bubble_pressure(methanol_water_mixture, thermo, opts):
    P, y, conv_flag = bubblePressure(methanol_water_mixture, thermo, opts)

    assert conv_flag == 1, f"bubblePressure did not converge"
    assert P > 0, f"bubble pressure negative: {P}"
    assert np.isfinite(P)
    assert np.all(np.isfinite(y))
    assert abs(np.sum(y) - 1.0) < 1e-6

    z = methanol_water_mixture.mole_fraction
    assert y[0] > z[0], f"y_CH4O ({y[0]}) should be > z_CH4O ({z[0]})"


def test_bubble_temperature(methanol_water_mixture, thermo, opts):
    T, y, conv_flag = bubbleTemperature(methanol_water_mixture, thermo, opts)

    assert conv_flag == 1, f"bubbleTemperature did not converge"
    assert T > 0, f"bubble temperature negative: {T}"
    assert np.isfinite(T)
    assert np.all(np.isfinite(y))
    assert abs(np.sum(y) - 1.0) < 1e-6


def test_dew_pressure(methanol_water_mixture, thermo, opts):
    P, x, conv_flag = dewPressure(methanol_water_mixture, thermo, opts)

    assert conv_flag == 1, f"dewPressure did not converge"
    assert P > 0, f"dew pressure negative: {P}"
    assert np.isfinite(P)
    assert np.all(np.isfinite(x))
    assert abs(np.sum(x) - 1.0) < 1e-6


def test_dew_temperature(methanol_water_mixture, thermo, opts):
    T, x, conv_flag = dewTemperature(methanol_water_mixture, thermo, opts)

    assert conv_flag == 1, f"dewTemperature did not converge"
    assert T > 0, f"dew temperature negative: {T}"
    assert np.isfinite(T)
    assert np.all(np.isfinite(x))
    assert abs(np.sum(x) - 1.0) < 1e-6
