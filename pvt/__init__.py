"""PVTtool — Thermodynamic phase-equilibrium calculations using cubic equations of state.

Provides VLE/LLE flash, phase stability testing, bubble/dew point calculations,
and residual property evaluation for multicomponent mixtures.

Usage:
    from pvt import Component, Mixture, ThermoModel, PREOS
    from pvt import vleflash, stabilityTest, FlashOptions

    comps, flag = Component.from_database_array(['CH4', 'C2H6'])
    mix = Mixture(comps, 300, 5e6)
    mix.mole_fraction = [0.6, 0.4]
    thermo = ThermoModel()
    y, x, V = vleflash(mix, thermo)
"""

from .component import Component
from .bip import BIP
from .mixture import Mixture
from .thermo_model import ThermoModel
from .flash_options import FlashOptions

from .eos import (
    ResidualProps,
    select_z_roots,
    PREOS,
    SRKEOS,
    PR78EOS,
)

from .mixing_rules import mixing_rule

from .activity import NRTL, Wilson, UNIQUAC, Margules2

from .utils import mynormalize, bisection

from .kvalues import kval_estimate, kvalue, kvalueLLE, fugacity

from .rachford_rice import RachfordRiceNR

from .flash import (
    massbalfunc,
    xy_calc,
    vleflash,
    vleflashnegative,
    lleflash,
    stabilityTest,
    stabilityLLETest,
)

from .saturation import (
    bubbleTemperature,
    bubblePressure,
    dewTemperature,
    dewPressure,
)

# Convenience wrappers matching MATLAB Tools/ functions


def add_components(names: list[str]) -> tuple[list[Component], list[int]]:
    """Load Component objects from the built-in database.

    Convenience wrapper around Component.from_database_array.

    Args:
        names: List of component names or formulas.

    Returns:
        Tuple of (components, not_found_indices).
    """
    return Component.from_database_array(names)


def add_mixture(
    components: list[Component], T_K: float, p_Pa: float
) -> Mixture:
    """Create a Mixture object from components, temperature, and pressure.

    Convenience wrapper around Mixture constructor.

    Args:
        components: List of Component objects.
        T_K: Temperature [K].
        p_Pa: Pressure [Pa].

    Returns:
        Mixture object with equimolar composition and zeroed BIP.
    """
    return Mixture(components, T_K, p_Pa)


def add_thermo() -> ThermoModel:
    """Create a ThermoModel with default settings.

    Convenience wrapper around ThermoModel constructor.

    Returns:
        ThermoModel with PR EOS, NRTL activity model, vdW mixing, liquid phase.
    """
    return ThermoModel()


def zero_bip(components: list[Component]) -> BIP:
    """Create a zero-initialised BIP for a given set of components.

    Convenience wrapper around BIP constructor.

    Args:
        components: List of Component objects.

    Returns:
        BIP with all matrices initialised to zero.
    """
    return BIP(len(components))
