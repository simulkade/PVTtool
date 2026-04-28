"""Activity coefficient models for PVTtool.

All models follow the signature:
    (temperature, x, components, bip) -> (gErt, gama)

where gErt = g^E/(R*T) and gama is a 1D array of activity coefficients.

Models: NRTL, Wilson, UNIQUAC, Margules2
"""

import numpy as np

from ._constants import R
from .bip import BIP
from .component import Component


def NRTL(
    temperature: float,
    x: np.ndarray,
    components: list[Component],
    bip: BIP,
) -> tuple[float, np.ndarray]:
    """Non-Random Two-Liquid (NRTL) activity coefficient model.

    g^E/(RT) = sum_i x_i * (sum_j tau_ji*G_ji*x_j) / (sum_k G_ki*x_k)

    tau_ij = A_ij/(R*T), G_ij = exp(-alfa_ij * tau_ij)
    A_ij = Acons + Atdep*T + Atdep2*T² + Atdepm1/T + Atdeplog*ln(T)

    Args:
        temperature: Temperature [K].
        x: Mole fractions (1D array).
        components: List of Component objects.
        bip: BIP object with NRTL fields.

    Returns:
        Tuple of (gErt, gama) where gErt = g^E/(R*T), gama is activity coefficients.
    """
    T = temperature
    N = len(x)

    # Build A_ij matrix
    A = (
        bip.NRTLcons
        + bip.NRTLtdep * T
        + bip.NRTLtdep2 * T**2
        + bip.NRTLtdepm1 / max(T, 1e-10)
        + bip.NRTLtdeplog * np.log(max(T, 1e-10))
    )
    tau = A / (R * T)
    G = np.exp(-bip.NRTLalfa * tau)

    if N == 2:
        # Explicit formulas for binary system
        a12 = tau[0, 1]
        a21 = tau[1, 0]
        b12 = G[0, 1]
        b21 = G[1, 0]

        gama_1 = np.exp(
            x[1] ** 2
            * (a21 * (b21 / (x[0] + x[1] * b21)) ** 2 + a12 * b12 / (x[1] + x[0] * b12) ** 2)
        )
        gama_2 = np.exp(
            x[0] ** 2
            * (a12 * (b12 / (x[1] + x[0] * b12)) ** 2 + a21 * b21 / (x[0] + x[1] * b21) ** 2)
        )
        gama = np.array([gama_1, gama_2])
        gErt = float(x[0] * np.log(gama_1) + x[1] * np.log(gama_2))
        return gErt, gama

    # General N-component case (vectorized)
    xG = x @ G  # 1xN
    # sum_j x_j * tau_ji * G_ji
    num = x @ (G * tau.T)
    gErt = float(np.dot(x, num / xG))

    # Activity coefficients
    part1 = num / xG
    part2_num = (x / xG) @ (G * tau).T
    part3_num = (x / xG**2 * (x @ (tau * G))) @ G.T
    gama = np.exp(part1 + part2_num - part3_num)

    return gErt, gama


def Wilson(
    temperature: float,
    x: np.ndarray,
    components: list[Component],
    bip: BIP,
) -> tuple[float, np.ndarray]:
    """Wilson activity coefficient model.

    A_ij = Wilsoncons + Wilsontdep * T
    Lambda_ij = Vj/Vi * exp(-A_ij/(R*T))

    Note: currently a stub — returns zeros.

    Args:
        temperature: Temperature [K].
        x: Mole fractions (1D array).
        components: List of Component objects.
        bip: BIP object with Wilson fields.

    Returns:
        Tuple of (gErt, gama).
    """
    N = len(x)
    gErt = 0.0
    gama = np.ones(N)
    return gErt, gama


def UNIQUAC(
    temperature: float,
    x: np.ndarray,
    components: list[Component],
    bip: BIP,
) -> tuple[float, np.ndarray]:
    """UNIQUAC (Universal Quasi-Chemical) activity coefficient model.

    Note: stub — returns ones.

    Args:
        temperature: Temperature [K].
        x: Mole fractions (1D array).
        components: List of Component objects.
        bip: BIP object with UNIQUAC fields.

    Returns:
        Tuple of (gErt, gama).
    """
    N = len(x)
    gErt = 0.0
    gama = np.ones(N)
    return gErt, gama


def Margules2(
    temperature: float,
    x: np.ndarray,
    components: list[Component],
    bip: BIP,
) -> tuple[float, np.ndarray]:
    """Two-parameter Margules activity coefficient model.

    Note: stub — returns ones.

    Args:
        temperature: Temperature [K].
        x: Mole fractions (1D array).
        components: List of Component objects.
        bip: BIP object.

    Returns:
        Tuple of (gErt, gama).
    """
    N = len(x)
    gErt = 0.0
    gama = np.ones(N)
    return gErt, gama
