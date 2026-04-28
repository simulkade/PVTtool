"""Correlation equation implementations for pure-component properties.

Each equation type corresponds to a string identifier from the database.
All functions accept (T, coefs) or (Tr, coefs) and return a scalar value in SI units.
"""

import numpy as np


def psat_antoine(T: float, coefs: np.ndarray) -> float:
    """Extended Antoine/DIPPR-101 vapor pressure equation.

    ln(P_sat) = c1 + c2/T + c3*ln(T) + c4*T^c5  [Pa]
    """
    c1, c2, c3, c4, c5 = coefs
    return float(np.exp(c1 + c2 / T + c3 * np.log(T) + c4 * T**c5))


def dhvap_dippr106(Tr: float, coefs: np.ndarray) -> float:
    """DIPPR-106 enthalpy of vaporization correlation.

    dh_vap = c1 * (1 - Tr)^(c2 + c3*Tr + c4*Tr^2) / 1000  [J/mol]
    Note: MATLAB version divides by 1000 (conversion from J/kmol to kJ/mol? check).
    Actually the MATLAB stores this in J/mol — coefficients are already scaled.

    We match the MATLAB behavior: the raw coefs give J/mol after this formula.
    """
    c1, c2, c3, c4 = coefs
    exponent = c2 + c3 * Tr + c4 * Tr**2
    return float(c1 * (1.0 - Tr) ** exponent)


def cpliq_poly(T: float, coefs: np.ndarray) -> float:
    """Polynomial liquid heat capacity equation.

    cp = (c1 + c2*T + c3*T^2 + c4*T^3 + c5*T^4) / 1000  [J/(mol·K)]
    The division by 1000 matches the MATLAB function handle convention.
    """
    c1, c2, c3, c4, c5 = coefs
    return float((c1 + c2 * T + c3 * T**2 + c4 * T**3 + c5 * T**4) / 1000.0)


def cpliq_dippr100(Tr: float, coefs: np.ndarray) -> float:
    """DIPPR-100 liquid heat capacity equation.

    cp = (c1^2/(1-Tr) + c2 - 2*c1*c3*(1-Tr) - c1*c4*(1-Tr)^2
         - c3^2/3*(1-Tr)^3 - c3*c4/2*(1-Tr)^4 - c4^2/5*(1-Tr)^5) / 1000

    Note: MATLAB expression has c3^(2/3) which should be c3^2/3.
    """
    c1, c2, c3, c4, c5 = coefs
    t = 1.0 - Tr
    val = (
        c1**2 / t
        + c2
        - 2.0 * c1 * c3 * t
        - c1 * c4 * t**2
        - c3**2 / 3.0 * t**3
        - c3 * c4 / 2.0 * t**4
        - c4**2 / 5.0 * t**5
    )
    return float(val / 1000.0)


def cpig_dippr107(T: float, coefs: np.ndarray) -> float:
    """Aly-Lee / DIPPR-107 ideal-gas heat capacity equation.

    cp = (c1 + c2*(c3/T / sinh(c3/T))^2 + c4*(c5/T / cosh(c5/T))^2) / 1000
    """
    c1, c2, c3, c4, c5 = coefs
    if T <= 0:
        return 0.0
    x3 = c3 / T
    x5 = c5 / T
    val = c1 + c2 * (x3 / np.sinh(x3)) ** 2 + c4 * (x5 / np.cosh(x5)) ** 2
    return float(val / 1000.0)


def cpig_poly(T: float, coefs: np.ndarray) -> float:
    """Polynomial ideal-gas heat capacity equation.

    cp = (c1 + c2*T + c3*T^2 + c4*T^3 + c5*T^4) / 1000
    """
    c1, c2, c3, c4, c5 = coefs
    return float((c1 + c2 * T + c3 * T**2 + c4 * T**3 + c5 * T**4) / 1000.0)


def cpig_loglin(T: float, coefs: np.ndarray) -> float:
    """Log-linear ideal-gas heat capacity equation.

    cp = (c1 + c2*ln(T) + c3/T + c4*T) / 1000
    """
    c1, c2, c3, c4, c5 = coefs
    return float((c1 + c2 * np.log(T) + c3 / T + c4 * T) / 1000.0)


# Map equation type strings to functions
EQUATION_FUNCTIONS = {
    "psat_antoine": psat_antoine,
    "dhvap_dippr106": dhvap_dippr106,
    "cpliq_poly": cpliq_poly,
    "cpliq_dippr100": cpliq_dippr100,
    "cpig_dippr107": cpig_dippr107,
    "cpig_poly": cpig_poly,
    "cpig_loglin": cpig_loglin,
}
