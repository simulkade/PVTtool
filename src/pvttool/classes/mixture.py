"""Multicomponent mixture state."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from .bip import BIP
from .component import Component


@dataclass
class Mixture:
    """Multicomponent mixture at fixed temperature and pressure.

    Parameters
    ----------
    components : list[Component]
        Pure-component objects in the mixture.
    temperature : float
        Temperature [K].
    pressure : float
        Pressure [Pa].

    Attributes
    ----------
    mole_fraction : np.ndarray, shape (n,)
        Overall mole fractions (initialised equimolar; must sum to 1).
    bip : BIP
        Binary interaction parameters (all zeros by default).
    """

    components: list[Component]
    temperature: float
    pressure: float
    mole_fraction: np.ndarray = field(init=False)
    bip: BIP = field(init=False)

    def __post_init__(self):
        n = len(self.components)
        self.mole_fraction = np.ones(n) / n
        self.bip = BIP(n)
