"""Tests for dataclass constructors and database loading."""

import numpy as np
import pytest

from pvttool import BIP, Component, FlashOptions, Mixture, ThermoModel, add_components, preos


def test_component_from_database_by_formula():
    comp = Component.from_database("CH4")
    assert comp.name.lower() in ("methane",)
    assert comp.formula == "CH4"
    assert comp.Tc == pytest.approx(190.564, rel=1e-3)
    assert comp.Pc == pytest.approx(4.59e6, rel=1e-2)
    assert comp.acentric_factor == pytest.approx(0.0115, abs=0.002)


def test_component_from_database_by_name():
    comp = Component.from_database("Methane")
    assert comp.formula == "CH4"


def test_component_from_database_case_insensitive():
    comp = Component.from_database("methane")
    assert comp.formula == "CH4"


def test_component_not_found():
    with pytest.raises(KeyError):
        Component.from_database("NotARealComponent_xyz")


def test_from_database_array():
    comps = Component.from_database_array(["Methane", "Ethane"])
    assert len(comps) == 2
    assert comps[0].formula == "CH4"
    assert comps[1].formula == "C2H6"


def test_component_vapor_pressure():
    comp = Component.from_database("Methane")
    psat = comp.vapor_pressure(150.0)
    assert psat > 0
    assert 1e4 < psat < 1e7  # physically plausible range for methane at 150 K


def test_component_cp_ig():
    comp = Component.from_database("Methane")
    cp = comp.cp_ig(300.0)
    assert cp > 0
    assert 20.0 < cp < 100.0  # J/(mol*K), methane ~35 J/(mol*K) at 300 K


def test_bip_initialization():
    bip = BIP(3)
    assert bip.eos_cons.shape == (3, 3)
    assert bip.eos_tdep.shape == (3, 3)
    assert bip.nrtl_cons.shape == (3, 3)
    assert bip.nrtl_alfa.shape == (3, 3)
    assert bip.uniquac_r.shape == (3,)
    assert np.all(bip.eos_cons == 0.0)
    assert np.all(bip.nrtl_cons == 0.0)


def test_mixture_equimolar_default():
    comps = add_components(["Methane", "Ethane"])
    mix = Mixture(comps, 300.0, 1e6)
    assert len(mix.components) == 2
    assert mix.mole_fraction.sum() == pytest.approx(1.0)
    assert mix.mole_fraction[0] == pytest.approx(0.5)
    assert mix.bip.eos_cons.shape == (2, 2)


def test_thermo_model_defaults():
    thermo = ThermoModel()
    assert thermo.eos is preos
    assert thermo.mixing_rule == 1
    assert thermo.phase == 1
    assert thermo.fugacity_switch == 1


def test_flash_options_defaults():
    opts = FlashOptions()
    assert opts.accuracy == pytest.approx(1e-7)
    assert opts.iteration == 100
    assert opts.trivial_solution_max_error == pytest.approx(1e-5)
    assert opts.convergence_max_error == pytest.approx(1e-10)
    assert opts.max_iteration == 50
