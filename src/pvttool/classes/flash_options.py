"""Flash calculation convergence options."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class FlashOptions:
    """Convergence tolerances and iteration limits for flash calculations.

    Attributes
    ----------
    accuracy : float
        Convergence tolerance for VLE/LLE flash (default 1e-7).
    iteration : int
        Maximum successive-substitution iterations for flash (default 100).
    trivial_solution_max_error : float
        Stability test: sum(log(K)^2) threshold for trivial solution
        detection (default 1e-5).
    convergence_max_error : float
        Stability test: sum((Ri-1)^2) threshold for convergence
        (default 1e-10).
    max_iteration : int
        Maximum iterations for stability test (default 50).
    """

    accuracy: float = 1e-7
    iteration: int = 100
    trivial_solution_max_error: float = 1e-5
    convergence_max_error: float = 1e-10
    max_iteration: int = 50
