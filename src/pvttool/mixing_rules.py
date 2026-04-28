"""Mixing rules for cubic equations of state."""

from __future__ import annotations

import numpy as np

from pvttool.classes.mixture import Mixture
from pvttool.classes.thermo_model import ThermoModel


_R = 8.314


def mixing_rule(
    mixture: Mixture,
    thermo: ThermoModel,
    ai: np.ndarray,
    bi: np.ndarray,
    s1: float,
    Q: np.ndarray,
) -> tuple[float, float]:
    """Compute mixture EOS parameters a and b using the specified mixing rule.

    Parameters
    ----------
    mixture : Mixture
        Mixture state (provides temperature, mole_fraction, bip, components).
    thermo : ThermoModel
        Thermodynamic model (provides mixing_rule number and activity_model).
    ai : np.ndarray, shape (n,)
        Pure-component a(T) parameters [Pa*m^6/mol^2].
    bi : np.ndarray, shape (n,)
        Pure-component b parameters [m^3/mol].
    s1 : float
        Huron-Vidal constant (0.623225 for PR, ln(2) for SRK).
    Q : np.ndarray, shape (2,)
        MHV constants [q1, q2]; zeros for vdW and Huron-Vidal.

    Returns
    -------
    a : float
        Mixture a parameter.
    b : float
        Mixture b parameter.

    Raises
    ------
    NotImplementedError
        If mixing_rule == 5 (Wong-Sandler, not implemented).
    """
    T = mixture.temperature
    x = mixture.mole_fraction
    bip = mixture.bip
    mr = thermo.mixing_rule
    N = len(x)

    b = float(np.dot(x, bi))

    if mr == 1:
        kij = bip.eos_cons + bip.eos_tdep * T
        a = 0.0
        for i in range(N):
            for j in range(N):
                a += x[i] * x[j] * np.sqrt(ai[i] * ai[j]) * (1.0 - kij[i, j])
        return a, b

    activity_fn = thermo.activity_model
    gErt, _ = activity_fn(T, x, mixture.components, bip)

    if mr == 2:
        a = b * (float(np.dot(x, ai / bi)) - gErt * _R * T / s1)
        return a, b

    q1, q2 = float(Q[0]), float(Q[1])
    alphai = ai / (bi * _R * T)

    if mr == 3:
        alpha = (gErt + float(np.dot(x, np.log(b / bi))) + q1 * float(np.dot(x, alphai))) / q1
        a = float(alpha) * b * _R * T
        return a, b

    if mr == 4:
        C = -(gErt + float(np.dot(x, np.log(b / bi)))
              + q1 * float(np.dot(x, alphai))
              + q2 * float(np.dot(x, alphai**2)))
        B_mhv = q1
        A_mhv = q2
        disc = B_mhv**2 - 4.0 * A_mhv * C
        alpha = (-B_mhv - np.sqrt(max(disc, 0.0))) / (2.0 * A_mhv)
        a = float(np.real(alpha)) * b * _R * T
        return a, b

    raise NotImplementedError("Wong-Sandler mixing rule (rule 5) is not yet implemented.")
