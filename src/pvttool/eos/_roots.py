"""Robust Z-root selection for cubic equations of state."""

from __future__ import annotations

import numpy as np


def select_z_roots(z_roots: np.ndarray, B_coef: float) -> tuple[float, float]:
    """Select liquid and vapor compressibility factors from cubic EOS roots.

    Filters the three roots returned by ``np.roots`` for physically meaningful
    values:

    * Roots whose imaginary-to-real ratio exceeds ``1e-6`` are rejected.
    * Roots with real part <= ``B_coef`` are rejected (required so that
      ``Z - B > 0`` and ``ln(Z - B)`` is defined).

    When only one physical root survives (supercritical or single-phase),
    ``z_liq == z_vap``.

    Parameters
    ----------
    z_roots : np.ndarray, shape (3,)
        Complex roots of the cubic EOS polynomial.
    B_coef : float
        Reduced co-volume ``B = b * P / (R * T)``.

    Returns
    -------
    z_liq : float
        Liquid compressibility factor (smallest physical root).
    z_vap : float
        Vapor compressibility factor (largest physical root).
    """
    tol = 1e-6
    real_z = np.real(z_roots)
    is_real = np.abs(np.imag(z_roots)) / np.maximum(np.abs(real_z), 1.0) < tol
    is_physical = is_real & (real_z > B_coef)
    phys_z = np.sort(real_z[is_physical])

    if len(phys_z) >= 2:
        return float(phys_z[0]), float(phys_z[-1])
    if len(phys_z) == 1:
        return float(phys_z[0]), float(phys_z[0])

    # Fallback: root with smallest imaginary part
    idx = int(np.argmin(np.abs(np.imag(z_roots))))
    zL = float(np.real(z_roots[idx]))
    if zL <= B_coef:
        zL = B_coef + 1e-8
    return zL, zL
