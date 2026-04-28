"""Mixing rules for cubic equations of state.

Computes mixture parameters a and b from pure-component parameters.
Supports van der Waals (1), Huron-Vidal (2), MHV1 (3), MHV2 (4).
Rule 5 (Wong-Sandler) is deliberately incomplete.
"""

import numpy as np

from ._constants import R
from .mixture import Mixture
from .thermo_model import ThermoModel


def mixing_rule(
    mixture: Mixture,
    thermo: ThermoModel,
    ai: np.ndarray,
    bi: np.ndarray,
    s1: float,
    Q: np.ndarray,
) -> tuple[float, float]:
    """Compute mixture a and b parameters for a cubic EOS.

    Args:
        mixture: Mixture object with components, mole_fraction, temperature, bip.
        thermo: ThermoModel with mixingrule, activity_model.
        ai: Pure-component a_i parameter array.
        bi: Pure-component b_i parameter array.
        s1: Mixing rule constant (0.623225 for PR, ln(2) for SRK).
        Q: MHV parameter array [q1, q2] if applicable.

    Returns:
        Tuple of (a, b) mixture parameters.

    Raises:
        NotImplementedError: If mixing rule 5 (Wong-Sandler) is requested.
    """
    rule = thermo.mixingrule
    x = mixture.mole_fraction
    T = mixture.temperature
    BIP = mixture.bip
    n = len(x)

    if rule == 1:
        # van der Waals one-fluid mixing rule
        b = float(np.dot(x, bi))
        a_matrix = np.sqrt(np.outer(ai, ai)) * (1.0 - BIP.EOScons - BIP.EOStdep * T)
        a = float(np.dot(x, np.dot(a_matrix, x)))
        return a, b

    elif rule in (2, 3, 4):
        # Activity coefficient based mixing rules (HV, MHV1, MHV2)
        # Requires activity model
        activity_fn = thermo.activity_model
        gErt, _ = activity_fn(T, x, mixture.components, BIP)
        b = float(np.dot(x, bi))

        if rule == 2:
            # Huron-Vidal
            a_hv = float(np.dot(x, ai / bi))
            a = float(b * (a_hv - gErt * R * T / s1))
            return a, b

        elif rule == 3:
            # MHV1
            q1 = Q[0]
            alpha_i = ai / (bi * R * T)
            alpha_avg = float(np.dot(x, alpha_i))
            sum_ln = float(np.dot(x, np.log(b / bi)))
            alpha = (gErt + sum_ln + q1 * alpha_avg) / q1
            a = float(alpha * b * R * T)
            return a, b

        elif rule == 4:
            # MHV2 — solve quadratic: q2*alpha² + q1*alpha + C = 0
            q1 = Q[0]
            q2 = Q[1]
            alpha_i = ai / (bi * R * T)
            alpha_avg = float(np.dot(x, alpha_i))
            alpha2_avg = float(np.dot(x, alpha_i**2))
            sum_ln = float(np.dot(x, np.log(b / bi)))
            C = -gErt - sum_ln - q1 * alpha_avg - q2 * alpha2_avg
            discriminant = q1**2 - 4.0 * q2 * C
            # Select the smaller (more negative) root, consistent with MATLAB
            alpha = (-q1 - np.sqrt(complex(discriminant))) / (2.0 * q2)
            alpha = float(np.real(alpha))
            a = float(alpha * b * R * T)
            return a, b

    elif rule == 5:
        # Wong-Sandler — deliberately incomplete
        raise NotImplementedError(
            "Wong-Sandler mixing rule (rule 5) is not yet implemented."
        )

    raise ValueError(f"Unknown mixing rule: {rule}")
