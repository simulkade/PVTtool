"""Unit tests for BIP, Component, Mixture, ThermoModel, FlashOptions.

Python equivalent of Tests/test_classes.m.
"""

import numpy as np
import pytest

from pvt import BIP, Component, Mixture, ThermoModel, FlashOptions
from pvt import add_components, add_mixture, add_thermo
from pvt.eos import PREOS


class TestBIP:
    def test_bip2_sizes(self):
        b2 = BIP(2)
        assert b2.EOScons.shape == (2, 2)
        assert b2.NRTLalfa.shape == (2, 2)
        assert b2.UNIQUACR.shape == (2,)  # 1D length 2

    def test_bip2_zero_init(self):
        b2 = BIP(2)
        assert np.all(b2.EOScons == 0)
        assert np.all(b2.NRTLcons == 0)

    def test_bip3_uniquac_cons(self):
        b3 = BIP(3)
        assert b3.UNIQUACcons.shape == (3, 3)


class TestComponentFromDatabaseArray:
    def test_ch4_c2h6(self):
        comps, flag = Component.from_database_array(["CH4", "C2H6"])
        assert flag == []
        assert len(comps) == 2
        assert all(isinstance(c, Component) for c in comps)
        assert len(comps[0].name) > 0
        assert comps[0].Tc > 0
        assert comps[0].Pc > 0
        assert comps[0].acentric_factor >= 0


class TestComponentFromDatabase:
    def test_h2o(self):
        c = Component.from_database("H2O")
        assert isinstance(c, Component)
        assert c.Tc > 0


class TestAddComponents:
    def test_ch4_h2o(self):
        comps, flag = add_components(["CH4", "H2O"])
        assert flag == []
        assert len(comps) == 2
        assert all(isinstance(c, Component) for c in comps)

    def test_missing_component(self):
        _, flag = add_components(["ThisDoesNotExist999"])
        assert len(flag) > 0


class TestMixture:
    def test_constructor(self):
        comps, _ = add_components(["CH4", "C2H6"])
        mix = Mixture(comps, 300, 1e6)
        assert isinstance(mix, Mixture)
        assert mix.temperature == 300
        assert mix.pressure == 1e6
        assert len(mix.mole_fraction) == 2
        assert abs(np.sum(mix.mole_fraction) - 1) < 1e-12
        assert isinstance(mix.bip, BIP)


class TestAddMixture:
    def test_wrapper(self):
        comps, _ = add_components(["CH4", "C2H6"])
        mix = add_mixture(comps, 350, 2e6)
        assert isinstance(mix, Mixture)
        assert mix.temperature == 350


class TestThermoModel:
    def test_defaults(self):
        th = ThermoModel()
        assert isinstance(th, ThermoModel)
        assert th.mixingrule == 1
        assert th.fugacity_switch == 1
        assert th.EOS is PREOS


class TestAddThermo:
    def test_wrapper(self):
        th = add_thermo()
        assert isinstance(th, ThermoModel)


class TestFlashOptions:
    def test_defaults(self):
        opts = FlashOptions()
        assert isinstance(opts, FlashOptions)
        assert opts.accuracy == 1e-7
        assert opts.iteration == 100
        assert opts.convergenceMaxError == 1e-10
