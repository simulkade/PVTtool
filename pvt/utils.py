"""Utility functions for PVTtool.

Provides normalisation and bisection root-finding helpers.
"""

import numpy as np
from collections.abc import Callable


def mynormalize(a: np.ndarray) -> np.ndarray:
    """Normalise a vector so its elements sum to 1.

    Args:
        a: 1D numpy array.

    Returns:
        Normalised 1D array a / sum(a).
    """
    s = np.sum(a)
    if s == 0:
        return np.zeros_like(a)
    return a / s


def bisection(f: Callable[[float], float], a: float, b: float, tol: float = 1e-6, max_iter: int = 100) -> float:
    """Find a root of f(x) = 0 in interval [a, b] using the bisection method.

    Args:
        f: Function of a single variable.
        a: Left bracket.
        b: Right bracket.
        tol: Convergence tolerance on |f(x)|.
        max_iter: Maximum number of iterations.

    Returns:
        Root approximation.

    Raises:
        RuntimeError: If f(a) and f(b) have the same sign (no root guaranteed).
        RuntimeError: If maximum iterations reached without convergence.
    """
    fa = f(a)
    fb = f(b)

    if fa * fb > 0:
        raise RuntimeError(
            f"bisection: f(a) and f(b) have same sign ({fa}, {fb}). "
            f"No root is guaranteed in [{a}, {b}]."
        )

    x_left = a
    x_right = b
    f_left = fa

    for _ in range(max_iter):
        s = (x_left + x_right) / 2.0
        fs = f(s)
        if abs(fs) < tol:
            return s
        if f_left * fs < 0:
            x_right = s
        else:
            x_left = s
            f_left = fs

    return (x_left + x_right) / 2.0
