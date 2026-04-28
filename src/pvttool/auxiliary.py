"""Auxiliary utilities."""

from __future__ import annotations

import numpy as np

from pvttool.classes.mixture import Mixture
from pvttool.flash._kvalue import kval_estimate


def wilson_correlation(mixture: Mixture) -> np.ndarray:
    """Wilson correlation K-value estimate (alias for kval_estimate).

    Parameters
    ----------
    mixture : Mixture
        Mixture state.

    Returns
    -------
    K : np.ndarray, shape (n,)
        Wilson K-values.
    """
    return kval_estimate(mixture)


def normalize(v: np.ndarray) -> np.ndarray:
    """Normalize a vector to sum to 1.

    Parameters
    ----------
    v : np.ndarray
        Input vector.

    Returns
    -------
    np.ndarray
        Normalized vector.
    """
    s = v.sum()
    return v / s if s > 0 else v
