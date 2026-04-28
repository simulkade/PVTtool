import numpy as np

from pvt import (
    Component,
    Mixture,
    PREOS,
    ResidualProps,
    SRKEOS,
    ThermoModel,
)

R = 8.314


def test_preos_ch4_liquid_props_negative():
    comps, _ = Component.from_database_array(["CH4"])
    mix = Mixture(comps, 150.0, 1e6)
    thermo = ThermoModel()
    thermo.phase = 1

    zL, zV, fug, HR, props = PREOS(mix, thermo, compute_props=True)

    assert isinstance(props, ResidualProps)
    assert props.SR < 0, f"Expected SR < 0, got {props.SR}"
    assert props.VR < 0, f"Expected VR < 0, got {props.VR}"
    assert props.HR < 0, f"Expected HR < 0, got {props.HR}"


def test_gr_identity():
    comps, _ = Component.from_database_array(["CH4"])
    mix = Mixture(comps, 150.0, 1e6)
    thermo = ThermoModel()
    thermo.phase = 1

    _, _, _, HR, props = PREOS(mix, thermo, compute_props=True)
    T = mix.temperature

    expected_GR = HR - T * props.SR
    assert abs(props.GR - expected_GR) < 1e-8, (
        f"GR identity failed: {props.GR} != {expected_GR}"
    )


def test_vr_identity():
    comps, _ = Component.from_database_array(["CH4"])
    mix = Mixture(comps, 150.0, 1e6)
    thermo = ThermoModel()
    thermo.phase = 1

    zL, zV, fug, HR, props = PREOS(mix, thermo, compute_props=True)
    T = mix.temperature
    P = mix.pressure

    expected_VR = R * T * (zL - 1.0) / P
    assert abs(props.VR - expected_VR) < 1e-8, (
        f"VR identity failed: {props.VR} != {expected_VR}"
    )


def test_srkeos_liquid_properties():
    comps, _ = Component.from_database_array(["CH4"])
    mix = Mixture(comps, 150.0, 1e6)
    thermo = ThermoModel()
    thermo.EOS = SRKEOS
    thermo.phase = 1

    zL, zV, fug, HR, props = SRKEOS(mix, thermo, compute_props=True)

    assert isinstance(props, ResidualProps)
    assert props.SR < 0, f"SRK: Expected SR < 0, got {props.SR}"
    assert props.VR < 0, f"SRK: Expected VR < 0, got {props.VR}"
    assert props.HR < 0, f"SRK: Expected HR < 0, got {props.HR}"

    T = mix.temperature
    assert abs(props.GR - (HR - T * props.SR)) < 1e-8
    assert abs(props.VR - R * T * (zL - 1.0) / mix.pressure) < 1e-8


def test_two_component_mixture_props():
    comps, _ = Component.from_database_array(["CH4", "C2H6"])
    mix = Mixture(comps, 200.0, 1e6)
    mix.mole_fraction = np.array([0.3, 0.7])
    thermo = ThermoModel()
    thermo.phase = 1

    zL, zV, fug, HR, props = PREOS(mix, thermo, compute_props=True)

    assert isinstance(props, ResidualProps)
    assert np.isfinite(props.SR)
    assert np.isfinite(props.VR)
    assert np.isfinite(props.HR)
    assert np.isfinite(props.GR)
    assert np.isfinite(props.Cp_R)
    assert np.isfinite(props.Cv_R)

    T = mix.temperature
    P = mix.pressure
    assert abs(props.GR - (HR - T * props.SR)) < 1e-8
    assert abs(props.VR - R * T * (zL - 1.0) / P) < 1e-8
