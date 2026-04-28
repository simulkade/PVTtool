"""Rachford-Rice equation solver and helpers."""

from __future__ import annotations

import numpy as np


def mass_bal_func(
    composition: np.ndarray,
    ki: np.ndarray,
    vapor_frac: float,
) -> tuple[float, float]:
    """Evaluate the Rachford-Rice equation and its derivative.

    f(V) = sum_i z_i*(K_i-1) / (1 + V*(K_i-1))
    df/dV = -sum_i z_i*(K_i-1)^2 / (1 + V*(K_i-1))^2

    Parameters
    ----------
    composition : np.ndarray, shape (n,)
        Overall mole fractions.
    ki : np.ndarray, shape (n,)
        K-values.
    vapor_frac : float
        Current vapor fraction estimate.

    Returns
    -------
    f : float
        Rachford-Rice residual.
    dfdv : float
        Derivative with respect to vapor fraction.
    """
    mask = composition != 0.0
    denom = 1.0 + vapor_frac * (ki[mask] - 1.0)
    km1 = ki[mask] - 1.0
    f = float(np.sum(composition[mask] * km1 / denom))
    dfdv = float(np.sum(-composition[mask] * km1**2 / denom**2))
    return f, dfdv


def xy_calc(
    composition: np.ndarray,
    vapor_frac: float,
    ki: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute liquid and vapor compositions from feed and K-values.

    Parameters
    ----------
    composition : np.ndarray, shape (n,)
        Overall mole fractions.
    vapor_frac : float
        Vapor mole fraction.
    ki : np.ndarray, shape (n,)
        K-values.

    Returns
    -------
    liquid_x : np.ndarray, shape (n,)
        Liquid mole fractions.
    vapor_y : np.ndarray, shape (n,)
        Vapor mole fractions.
    """
    liquid_x = composition / (1.0 + vapor_frac * (ki - 1.0))
    vapor_y = ki * liquid_x
    return liquid_x, vapor_y


def rachford_rice_nr(
    composition: np.ndarray,
    ki: np.ndarray,
    vapor_frac_est: float = 0.5,
) -> tuple[float, bool]:
    """Solve the Rachford-Rice equation by Newton-Raphson.

    Enforces physical bounds based on K-values and uses a relative
    convergence criterion.

    Parameters
    ----------
    composition : np.ndarray, shape (n,)
        Overall mole fractions.
    ki : np.ndarray, shape (n,)
        K-values.
    vapor_frac_est : float, optional
        Initial vapor fraction guess (default 0.5).

    Returns
    -------
    vf : float
        Converged vapor fraction.
    converged : bool
        True if Newton-Raphson converged within 100 iterations.
    """
    K_max = float(np.max(ki))
    K_min = float(np.min(ki))
    V_min = 1.0 / (1.0 - K_max) if K_max > 1.0 else -1e6
    V_max = 1.0 / (1.0 - K_min) if K_min < 1.0 else 1e6

    vf = float(np.clip(vapor_frac_est, V_min + 1e-10, V_max - 1e-10))
    converged = False
    for _ in range(100):
        f, dfdv = mass_bal_func(composition, ki, vf)
        if abs(dfdv) < 1e-300:
            break
        dvf = f / dfdv
        vf -= dvf
        vf = float(np.clip(vf, V_min + 1e-10, V_max - 1e-10))
        if abs(dvf) < 1e-10 * (1.0 + abs(vf)):
            converged = True
            break
    return vf, converged
