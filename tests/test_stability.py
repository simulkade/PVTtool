import numpy as np

from pvt import (
    Component,
    Mixture,
    ThermoModel,
    stabilityTest,
    stabilityLLETest,
)


def test_ch4_stable():
    comps, _ = Component.from_database_array(["CH4"])
    mix = Mixture(comps, 300.0, 1e5)
    thermo = ThermoModel()

    flag, SL, SV, result = stabilityTest(mix, thermo)

    assert result["overall"] == "stable", f"Expected stable, got {result['overall']}: {result['message']}"


def test_ch4_c10h22_unstable():
    comps, _ = Component.from_database_array(["CH4", "C10H22"])
    mix = Mixture(comps, 300.0, 5e6)
    mix.mole_fraction = np.array([0.5, 0.5])
    thermo = ThermoModel()

    flag, SL, SV, result = stabilityTest(mix, thermo)

    assert result["overall"] == "unstable", f"Expected unstable, got {result['overall']}: {result['message']}"


def test_h2o_c10h22_lle():
    comps, _ = Component.from_database_array(["H2O", "C10H22"])
    mix = Mixture(comps, 300.0, 100e5)
    mix.mole_fraction = np.array([0.5, 0.5])
    thermo = ThermoModel()

    flag, SL, SV, result = stabilityLLETest(mix, thermo)

    assert result["overall"] == "unstable", f"Expected LLE unstable, got {result['overall']}: {result['message']}"
