"""Rachford-Rice equation solver using bounded Newton-Raphson."""

import numpy as np


def RachfordRiceNR(
    composition: np.ndarray,
    ki: np.ndarray,
    vapor_frac_est: float,
) -> tuple[float, int]:
    """Solve the Rachford-Rice equation for vapor fraction using Newton-Raphson.

    Solves f(V) = sum_i z_i*(K_i-1) / (1 + V*(K_i-1)) = 0
    with bounds enforcement to keep V in the physical interval (Vmin, Vmax).

    Args:
        composition: [1 x N] feed mole fractions z_i.
        ki: [1 x N] K-values K_i.
        vapor_frac_est: Initial guess for vapor fraction V.

    Returns:
        Tuple of (vf, flag) where vf is the solved vapor fraction
        and flag is 1 if converged, 0 if max iterations reached.
    """
    Kmax = float(np.max(ki))
    Kmin = float(np.min(ki))

    # Physical bounds
    if Kmax > 1.0:
        Vmin = 1.0 / (1.0 - Kmax)
    else:
        Vmin = -1e10
    if Kmin < 1.0:
        Vmax = 1.0 / (1.0 - Kmin)
    else:
        Vmax = 1e10

    vapor_frac = min(max(vapor_frac_est, Vmin + 1e-10), Vmax - 1e-10)
    RRflag = 0

    for _ in range(100):
        f = 0.0
        dfdv = 0.0
        for i in range(len(composition)):
            if composition[i] != 0.0:
                denom = 1.0 + vapor_frac * (ki[i] - 1.0)
                f += composition[i] * (ki[i] - 1.0) / denom
                dfdv -= composition[i] * (ki[i] - 1.0) ** 2 / denom ** 2

        if abs(dfdv) < 1e-30:
            break

        dvf = f / dfdv
        vapor_frac = vapor_frac - dvf
        vapor_frac = min(max(vapor_frac, Vmin + 1e-10), Vmax - 1e-10)

        if abs(dvf) < 1e-10 * (1.0 + abs(vapor_frac)):
            RRflag = 1
            break

    return float(vapor_frac), RRflag
