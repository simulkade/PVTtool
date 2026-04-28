"""Wilson activity coefficient model."""

from __future__ import annotations

import numpy as np


def wilson(
    T: float,
    x: np.ndarray,
    components: list,
    bip,
) -> tuple[float, np.ndarray]:
    """Wilson activity coefficient model.

    Parameters
    ----------
    T : float
        Temperature [K].
    x : np.ndarray, shape (n,)
        Mole fractions.
    components : list[Component]
        Pure-component objects (not used; included for API uniformity).
    bip : BIP
        Binary interaction parameters.  Uses ``wilson_cons`` and
        ``wilson_tdep``.

    Returns
    -------
    g_ert : float
        Excess Gibbs energy divided by RT [dimensionless].
    gamma : np.ndarray, shape (n,)
        Activity coefficients.
    """
    A = bip.wilson_cons + bip.wilson_tdep * T
    N = len(x)

    if N == 2:
        x1, x2 = x[0], x[1]
        A12, A21 = A[0, 1], A[1, 0]
        g_ert = -x1 * np.log(x1 + A12 * x2) - x2 * np.log(x2 + A21 * x1)
        gamma = np.array([
            np.exp(-np.log(x1 + A12 * x2)
                   + x2 * (A12 / (x1 + A12 * x2) - A21 / (A21 * x1 + x2))),
            np.exp(-np.log(x2 + A21 * x1)
                   - x1 * (A12 / (x1 + A12 * x2) - A21 / (A21 * x1 + x2))),
        ])
        return float(g_ert), gamma

    xA = x @ A                       # shape (n,)  = sum_j x_j A_ji for each i... wait
    # Actually MATLAB: x*A' = (row vector)*(matrix) which sums over j for each i
    xAt = x @ A                      # x @ A  where A[i,j] is Lambda_ij
    g_ert = float(x @ np.log(xAt))
    part1 = -np.log(xAt) + 1.0
    part2 = -(x / xAt) @ A
    gamma = np.exp(part1 + part2)
    return g_ert, gamma
