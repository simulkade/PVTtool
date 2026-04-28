"""Convenience wrappers mirroring the MATLAB Tools/ functions."""

from __future__ import annotations

from pvttool.classes.bip import BIP
from pvttool.classes.component import Component
from pvttool.classes.flash_options import FlashOptions
from pvttool.classes.mixture import Mixture
from pvttool.classes.thermo_model import ThermoModel


def add_components(names: list[str]) -> list[Component]:
    """Load components from the built-in database.

    Parameters
    ----------
    names : list[str]
        Component names or chemical formulae (case-insensitive).

    Returns
    -------
    list[Component]
        Loaded components in the order of `names`.

    Raises
    ------
    KeyError
        If any name is not found in the database.
    """
    return Component.from_database_array(names)


def add_mixture(
    components: list[Component],
    T: float,
    p: float,
) -> Mixture:
    """Create a Mixture object with equimolar composition and zero BIP.

    Parameters
    ----------
    components : list[Component]
        Pure-component objects.
    T : float
        Temperature [K].
    p : float
        Pressure [Pa].

    Returns
    -------
    Mixture
        New mixture with equimolar ``mole_fraction`` and all-zero BIP.
    """
    return Mixture(components=components, temperature=T, pressure=p)


def add_thermo() -> ThermoModel:
    """Create a ThermoModel with default settings (PREOS, vdW mixing rule).

    Returns
    -------
    ThermoModel
        Default thermodynamic model configuration.
    """
    return ThermoModel()


def zero_bip(components: list[Component]) -> BIP:
    """Create a BIP object with all matrices initialised to zero.

    Parameters
    ----------
    components : list[Component]
        Components in the mixture (determines size n).

    Returns
    -------
    BIP
        All-zero binary interaction parameter container.
    """
    return BIP(n=len(components))
