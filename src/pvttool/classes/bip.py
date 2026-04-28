"""Binary interaction parameter container."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class BIP:
    """Binary interaction parameter matrices for an n-component mixture.

    All matrices are initialised to zeros.  Users set individual entries
    before calling flash or EOS routines.

    Parameters
    ----------
    n : int
        Number of components.

    Attributes
    ----------
    eos_cons : np.ndarray, shape (n, n)
        Constant part of EOS binary interaction parameter kij.
    eos_tdep : np.ndarray, shape (n, n)
        Temperature-dependent part of kij (kij = eos_cons + eos_tdep * T).
    nrtl_cons : np.ndarray, shape (n, n)
        NRTL constant interaction parameter A [J/mol].
    nrtl_tdep : np.ndarray, shape (n, n)
        NRTL linear T-coefficient B [J/(mol*K)].
    nrtl_tdep2 : np.ndarray, shape (n, n)
        NRTL quadratic T-coefficient C [J/(mol*K^2)].
    nrtl_tdepm1 : np.ndarray, shape (n, n)
        NRTL reciprocal T-coefficient D [J*K/mol].
    nrtl_tdeplog : np.ndarray, shape (n, n)
        NRTL log-T coefficient E [J/mol].
    nrtl_alfa : np.ndarray, shape (n, n)
        NRTL non-randomness parameter alpha [-].
    wilson_cons : np.ndarray, shape (n, n)
        Wilson constant interaction parameter [J/mol].
    wilson_tdep : np.ndarray, shape (n, n)
        Wilson temperature-dependent coefficient [J/(mol*K)].
    uniquac_cons : np.ndarray, shape (n, n)
        UNIQUAC constant interaction parameter [J/mol].
    uniquac_tdep : np.ndarray, shape (n, n)
        UNIQUAC temperature-dependent parameter.
    uniquac_r : np.ndarray, shape (n,)
        UNIQUAC volume (r) parameters per component.
    uniquac_q : np.ndarray, shape (n,)
        UNIQUAC surface (q) parameters per component.
    """

    n: int
    eos_cons: np.ndarray = field(init=False)
    eos_tdep: np.ndarray = field(init=False)
    nrtl_cons: np.ndarray = field(init=False)
    nrtl_tdep: np.ndarray = field(init=False)
    nrtl_tdep2: np.ndarray = field(init=False)
    nrtl_tdepm1: np.ndarray = field(init=False)
    nrtl_tdeplog: np.ndarray = field(init=False)
    nrtl_alfa: np.ndarray = field(init=False)
    wilson_cons: np.ndarray = field(init=False)
    wilson_tdep: np.ndarray = field(init=False)
    uniquac_cons: np.ndarray = field(init=False)
    uniquac_tdep: np.ndarray = field(init=False)
    uniquac_r: np.ndarray = field(init=False)
    uniquac_q: np.ndarray = field(init=False)

    def __post_init__(self):
        n = self.n
        self.eos_cons = np.zeros((n, n))
        self.eos_tdep = np.zeros((n, n))
        self.nrtl_cons = np.zeros((n, n))
        self.nrtl_tdep = np.zeros((n, n))
        self.nrtl_tdep2 = np.zeros((n, n))
        self.nrtl_tdepm1 = np.zeros((n, n))
        self.nrtl_tdeplog = np.zeros((n, n))
        self.nrtl_alfa = np.zeros((n, n))
        self.wilson_cons = np.zeros((n, n))
        self.wilson_tdep = np.zeros((n, n))
        self.uniquac_cons = np.zeros((n, n))
        self.uniquac_tdep = np.zeros((n, n))
        self.uniquac_r = np.zeros(n)
        self.uniquac_q = np.zeros(n)
