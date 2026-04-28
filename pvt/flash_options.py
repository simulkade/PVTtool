"""FlashOptions — convergence settings for flash and stability calculations."""

from dataclasses import dataclass


@dataclass
class FlashOptions:
    """Convergence settings for flash and stability calculations.

    All properties have sensible defaults; override as needed.

    Attributes:
        accuracy: Convergence tolerance for VLE/LLE flash (default 1e-7).
        iteration: Maximum successive substitution iterations (default 100).
        trivialSolutionMaxError: Threshold below which stability solution is trivial (default 1e-5).
        convergenceMaxError: Tolerance for convergence of stability loop (default 1e-10).
        maxIteration: Maximum stability test iterations (default 50).
    """

    accuracy: float = 1e-7
    iteration: int = 100
    trivialSolutionMaxError: float = 1e-5
    convergenceMaxError: float = 1e-10
    maxIteration: int = 50
