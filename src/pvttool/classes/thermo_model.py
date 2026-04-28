"""Thermodynamic model configuration."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional


@dataclass
class ThermoModel:
    """Configuration for the equation of state and mixing rule.

    Attributes
    ----------
    eos : Callable
        EOS function (default: preos).  Signature:
        ``(mixture, thermo) -> (z_liq, z_vap, fugacity_coefs, HR, props)``.
    activity_model : Callable
        Activity coefficient model (default: nrtl).  Signature:
        ``(T, x, components, bip) -> (g_ert, gamma)``.
    mixing_rule : int
        1 = van der Waals, 2 = Huron-Vidal, 3 = MHV1, 4 = MHV2.
    phase : int
        Which Z-root to use: 1 = liquid (smallest), 2 = vapor (largest).
    fugacity_switch : int
        1 = compute fugacity coefficients, 0 = skip (faster, returns zeros).
    """

    eos: Callable = field(default=None)
    activity_model: Callable = field(default=None)
    mixing_rule: int = 1
    phase: int = 1
    fugacity_switch: int = 1

    def __post_init__(self):
        # Defer import to avoid circular dependencies at module load time.
        if self.eos is None:
            from pvttool.eos.pr import preos
            self.eos = preos
        if self.activity_model is None:
            from pvttool.activity.nrtl import nrtl
            self.activity_model = nrtl
