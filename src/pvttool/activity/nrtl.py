"""Non-Random Two-Liquid (NRTL) activity coefficient model."""

from __future__ import annotations

import numpy as np


_R = 8.314


def nrtl(
    T: float,
    x: np.ndarray,
    components: list,
    bip,
) -> tuple[float, np.ndarray]:
    """Non-Random Two-Liquid activity coefficient model.

    Parameters
    ----------
    T : float
        Temperature [K].
    x : np.ndarray, shape (n,)
        Mole fractions.
    components : list[Component]
        Pure-component objects (not used by NRTL; included for API uniformity).
    bip : BIP
        Binary interaction parameters.  Relevant fields:
        ``nrtl_cons``, ``nrtl_tdep``, ``nrtl_tdep2``, ``nrtl_tdepm1``,
        ``nrtl_tdeplog``, ``nrtl_alfa``.

    Returns
    -------
    g_ert : float
        Excess Gibbs energy divided by RT [dimensionless].
    gamma : np.ndarray, shape (n,)
        Activity coefficients.
    """
    A = (bip.nrtl_cons
         + bip.nrtl_tdep * T
         + bip.nrtl_tdep2 * T**2
         + bip.nrtl_tdepm1 / T
         + bip.nrtl_tdeplog * np.log(T))
    alfa = bip.nrtl_alfa
    N = len(x)

    # Zero diagonal
    np.fill_diagonal(A, 0.0)

    taw = (A / (_R * T))
    G = np.exp(-alfa * taw)
    np.fill_diagonal(G, 0.0)

    if N == 2:
        t21, t12 = taw[1, 0], taw[0, 1]
        G21, G12 = G[1, 0], G[0, 1]
        x1, x2 = x[0], x[1]
        g_ert = x1 * x2 * (t21 * G21 / (x1 + x2 * G21) + t12 * G12 / (x2 + x1 * G12))
        gamma = np.array([
            np.exp(x2**2 * (t21 * (G21 / (x1 + x2 * G21))**2
                            + t12 * G12 / (x2 + x1 * G12)**2)),
            np.exp(x1**2 * (t12 * (G12 / (x2 + x1 * G12))**2
                            + t21 * G21 / (x1 + x2 * G21)**2)),
        ])
        return float(g_ert), gamma

    # General N-component case
    xG = x @ G                           # shape (n,)
    xGtaw = x @ (G * taw)               # shape (n,)
    g_ert = float(x @ (xGtaw / xG))

    part1 = xGtaw / xG
    part2 = (x / xG) @ (G * taw)
    part3 = (x / xG**2) * (x @ (taw * G)) @ G.T
    gamma = np.exp(part1 + part2 - part3)
    return g_ert, gamma
