"""BIP — Binary Interaction Parameters container.

Stores all matrices for binary interaction parameters used by equations of state
and activity coefficient models. All fields are zero-initialised for an
n-component mixture.

Matrices are n×n symmetric numpy arrays. UNIQUAC R and Q are 1D arrays of length n.
"""

from dataclasses import dataclass, field

import numpy as np


@dataclass
class BIP:
    """Binary Interaction Parameters for an n-component mixture.

    All interaction matrices are initialised to zero. Set non-zero values
    by assigning to individual fields after construction.

    EOS kij(T) = EOScons + EOStdep * T
    NRTL A(T) = NRTLcons + NRTLtdep*T + NRTLtdep2*T² + NRTLtdepm1/T + NRTLtdeplog*ln(T)
    Wilson A(T) = Wilsoncons + Wilsontdep*T
    UNIQUAC A(T) = UNIQUACcons + UNIQUACtdep*T

    Attributes:
        EOScons: Constant part of EOS kij.
        EOStdep: Temperature-dependent part of EOS kij.
        NRTLcons: Constant part of NRTL A_ij.
        NRTLtdep: B coefficient for NRTL (T term).
        NRTLtdep2: C coefficient for NRTL (T² term).
        NRTLtdepm1: D coefficient for NRTL (1/T term).
        NRTLtdeplog: E coefficient for NRTL (ln T term).
        NRTLalfa: Non-randomness parameter α_ij.
        Wilsoncons: Constant part of Wilson A_ij.
        Wilsontdep: Temperature-dependent part of Wilson A_ij.
        UNIQUACcons: Constant part of UNIQUAC A_ij.
        UNIQUACtdep: Temperature-dependent part of UNIQUAC A_ij.
        UNIQUACR: UNIQUAC volume parameter r_i (1D array).
        UNIQUACQ: UNIQUAC surface area parameter q_i (1D array).
    """

    EOScons: np.ndarray
    EOStdep: np.ndarray
    NRTLcons: np.ndarray
    NRTLtdep: np.ndarray
    NRTLtdep2: np.ndarray
    NRTLtdepm1: np.ndarray
    NRTLtdeplog: np.ndarray
    NRTLalfa: np.ndarray
    Wilsoncons: np.ndarray
    Wilsontdep: np.ndarray
    UNIQUACcons: np.ndarray
    UNIQUACtdep: np.ndarray
    UNIQUACR: np.ndarray
    UNIQUACQ: np.ndarray

    def __init__(self, n: int):
        """Initialise all BIP matrices to zero for n components.

        Args:
            n: Number of components in the mixture.
        """
        self.EOScons = np.zeros((n, n))
        self.EOStdep = np.zeros((n, n))
        self.NRTLcons = np.zeros((n, n))
        self.NRTLtdep = np.zeros((n, n))
        self.NRTLtdep2 = np.zeros((n, n))
        self.NRTLtdepm1 = np.zeros((n, n))
        self.NRTLtdeplog = np.zeros((n, n))
        self.NRTLalfa = np.zeros((n, n))
        self.Wilsoncons = np.zeros((n, n))
        self.Wilsontdep = np.zeros((n, n))
        self.UNIQUACcons = np.zeros((n, n))
        self.UNIQUACtdep = np.zeros((n, n))
        self.UNIQUACR = np.zeros(n)
        self.UNIQUACQ = np.zeros(n)
