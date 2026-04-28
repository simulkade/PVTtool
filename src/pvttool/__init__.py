"""pvttool — PVT flash calculation library.

Mirrors the MATLAB PVTtool API with Python snake_case naming.

Quick-start example
-------------------
>>> import numpy as np
>>> from pvttool import add_components, add_mixture, add_thermo, vle_flash
>>> comps = add_components(["Methanol", "Water"])
>>> mix = add_mixture(comps, T=350.0, p=101325.0)
>>> mix.mole_fraction = np.array([0.5, 0.5])
>>> y, x, V = vle_flash(mix, add_thermo())
"""

from pvttool.classes.component import Component
from pvttool.classes.bip import BIP
from pvttool.classes.mixture import Mixture
from pvttool.classes.thermo_model import ThermoModel
from pvttool.classes.flash_options import FlashOptions

from pvttool.eos.pr import preos
from pvttool.eos.srk import srkeos
from pvttool.eos.pr78 import pr78eos

from pvttool.activity.nrtl import nrtl
from pvttool.activity.wilson import wilson
from pvttool.activity.uniquac import uniquac

from pvttool.flash.vle import vle_flash, vle_flash_negative
from pvttool.flash.lle import lle_flash
from pvttool.flash.saturation import (
    bubble_pressure,
    bubble_temperature,
    dew_pressure,
    dew_temperature,
)
from pvttool.flash.stability import stability_test, stability_lle_test

from pvttool.tools import add_components, add_mixture, add_thermo, zero_bip
from pvttool.auxiliary import wilson_correlation, normalize

__all__ = [
    # Classes
    "Component", "BIP", "Mixture", "ThermoModel", "FlashOptions",
    # EOS
    "preos", "srkeos", "pr78eos",
    # Activity models
    "nrtl", "wilson", "uniquac",
    # Flash
    "vle_flash", "vle_flash_negative", "lle_flash",
    "bubble_pressure", "bubble_temperature",
    "dew_pressure", "dew_temperature",
    "stability_test", "stability_lle_test",
    # Tools
    "add_components", "add_mixture", "add_thermo", "zero_bip",
    # Utilities
    "wilson_correlation", "normalize",
]
