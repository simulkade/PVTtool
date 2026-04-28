"""ThermoModel — thermodynamic model configuration.

Specifies the equation of state, activity model, mixing rule, phase selection,
and whether to compute fugacity coefficients.
"""

from dataclasses import dataclass, field
from collections.abc import Callable


@dataclass
class ThermoModel:
    """Thermodynamic model configuration for flash and stability calculations.

    All fields have sensible defaults; override as needed.

    Attributes:
        EOS: Equation of state function (default: PREOS).
        activity_model: Activity coefficient model function (default: NRTL).
        mixingrule: Mixing rule 1=vdW, 2=HV, 3=MHV1, 4=MHV2.
        phase: Phase identifier 1=liquid (use Z_liq), 2=vapor (use Z_vap).
        fugacity_switch: 1 = compute fugacity coefficients, 0 = skip.
    """

    EOS: Callable = field(default=None)  # type: ignore[assignment]
    activity_model: Callable = field(default=None)  # type: ignore[assignment]
    mixingrule: int = 1
    phase: int = 1
    fugacity_switch: int = 1

    def __post_init__(self):
        """Set default function handles if not provided."""
        if self.EOS is None:
            from .eos import PREOS

            object.__setattr__(self, "EOS", PREOS)
        if self.activity_model is None:
            from .activity import NRTL

            object.__setattr__(self, "activity_model", NRTL)
