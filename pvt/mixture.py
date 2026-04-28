"""Mixture — multicomponent mixture at given temperature and pressure."""

from dataclasses import dataclass, field

import numpy as np

from .bip import BIP
from .component import Component


@dataclass
class Mixture:
    """Multicomponent mixture at specified temperature and pressure.

    Creates an equimolar mixture with zeroed BIP. Override mole_fraction
    and set bip fields after construction.

    Attributes:
        components: List of Component objects (1D).
        mole_fraction: Mole fractions (1D numpy array, sum = 1).
        pressure: System pressure [Pa].
        temperature: System temperature [K].
        bip: Binary Interaction Parameters object.
    """

    components: list[Component]
    temperature: float
    pressure: float
    mole_fraction: np.ndarray = field(default=None)  # type: ignore[assignment]
    bip: BIP = field(default=None)  # type: ignore[assignment]

    def __post_init__(self):
        """Initialise equimolar composition and zeroed BIP."""
        n = len(self.components)
        if self.mole_fraction is None:
            self.mole_fraction = np.ones(n) / n
        if self.bip is None:
            self.bip = BIP(n)
