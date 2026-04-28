"""Cubic equations of state for PVTtool.

Implements Peng-Robinson (1976), Soave-Redlich-Kwong (1972), and
Peng-Robinson with 1978 alpha correction. All functions return
compressibility factors, fugacity coefficients, residual enthalpy,
and optionally a full set of residual properties.

Includes the robust Z-root selection helper select_z_roots.
"""

from typing import NamedTuple

import numpy as np

from ._constants import R
from .mixing_rules import mixing_rule
from .mixture import Mixture
from .thermo_model import ThermoModel


class ResidualProps(NamedTuple):
    """Residual (departure) thermodynamic properties.

    Attributes:
        HR: Residual enthalpy [J/mol].
        SR: Residual entropy [J/(mol·K)].
        GR: Residual Gibbs free energy [J/mol].
        VR: Residual molar volume [m³/mol].
        Cp_R: Residual isobaric heat capacity [J/(mol·K)].
        Cv_R: Residual isochoric heat capacity [J/(mol·K)].
    """

    HR: float
    SR: float
    GR: float
    VR: float
    Cp_R: float
    Cv_R: float


def select_z_roots(z_roots: np.ndarray, B_coef: float) -> tuple[float, float]:
    """Select liquid and vapor compressibility factor roots from a cubic EOS.

    Filters out non-physical roots (complex, or real part <= B).

    Args:
        z_roots: Array of 3 complex roots from the cubic EOS.
        B_coef: Dimensionless B = bP/(RT).

    Returns:
        Tuple of (zL, zV) — smallest and largest physical real roots.
        If only one physical root, zL = zV.
    """
    tol = 1e-6
    real_z = np.real(z_roots)
    imag_z = np.abs(np.imag(z_roots))

    is_real = imag_z / np.maximum(np.abs(real_z), 1.0) < tol
    physical = is_real & (real_z > B_coef)
    phys_z = np.sort(real_z[physical])

    if len(phys_z) >= 2:
        return phys_z[0], phys_z[-1]
    elif len(phys_z) == 1:
        return phys_z[0], phys_z[0]
    else:
        idx = np.argmin(imag_z)
        z_fallback = real_z[idx]
        z_fallback = max(z_fallback, B_coef + 1e-10)
        return z_fallback, z_fallback

def PREOS(
    mixture: Mixture,
    thermo: ThermoModel,
    compute_props: bool = False,
) -> tuple[float, float, np.ndarray, float, ResidualProps | None]:
    """Peng-Robinson (1976) Equation of State.

    Supports van der Waals (rule 1), Huron-Vidal (rule 2), MHV1 (rule 3),
    and MHV2 (rule 4) mixing rules.

    Args:
        mixture: Mixture object with temperature, pressure, mole_fraction,
                 components, bip.
        thermo: ThermoModel with mixingrule, activity_model, phase,
                fugacity_switch.
        compute_props: If True, compute and return full residual properties.

    Returns:
        Tuple of (liquid_z, vapor_z, fugacity, HR, props).
        liquid_z: Liquid compressibility factor.
        vapor_z: Vapor compressibility factor.
        fugacity: [1 x N] fugacity coefficients (zero array if fugacity_switch == 0).
        HR: Residual molar enthalpy [J/mol].
        props: ResidualProps namedtuple if compute_props else None.
    """
    rule = thermo.mixingrule
    activity_fn = thermo.activity_model
    phase1 = thermo.phase
    fug_need = thermo.fugacity_switch

    critical_pres = np.array([c.Pc for c in mixture.components])
    critical_temp = np.array([c.Tc for c in mixture.components])
    omega = np.array([c.acentric_factor for c in mixture.components])
    BIP_data = mixture.bip
    x = mixture.mole_fraction
    p = mixture.pressure
    T = mixture.temperature
    N = len(critical_temp)

    fugacity = np.zeros(N)

    s1 = 0.623225  # Huron-Vidal constant for PR

    bi = 0.077796 * R * critical_temp / critical_pres
    aci = 0.457235 * (R * critical_temp) ** 2 / critical_pres
    mi = 0.37646 + (1.54226 - 0.26992 * omega) * omega
    Tr = T / critical_temp
    alfai = 1.0 + mi * (1.0 - np.sqrt(Tr))
    alfa = alfai ** 2
    ai = aci * alfa

    Q = np.array([0.0, 0.0])
    if rule == 3:
        Q = np.array([-0.53, 0.0])
    elif rule == 4:
        Q = np.array([-0.4347, -0.003654])

    a, b = mixing_rule(mixture, thermo, ai, bi, s1, Q)

    A_coef = a * p / (R * T) ** 2
    B_coef = b * p / (R * T)

    # PR cubic: Z^3 - (1-B)Z^2 + (A-3B^2-2B)Z - (AB-B^2-B^3) = 0
    poly_coef = np.array([
        1.0,
        -1.0 + B_coef,
        A_coef - B_coef * (2.0 + 3.0 * B_coef),
        -B_coef * (A_coef - B_coef * (1.0 + B_coef)),
    ])

    z_roots = np.roots(poly_coef)
    liquid_z, vapor_z = select_z_roots(z_roots, B_coef)

    if phase1 == 1:
        zz = liquid_z
    else:
        zz = vapor_z

    if fug_need == 1:
        TwoSqrt2 = 2.0 * np.sqrt(2.0)

        if rule == 1:
            part1 = bi / b * (zz - 1.0) - np.log(zz - B_coef)
            a_matrix = np.sqrt(np.outer(ai, ai)) * (
                1.0 - BIP_data.EOScons - BIP_data.EOStdep * T
            )
            part2 = x @ a_matrix
            part3 = (
                A_coef
                / (2.828 * B_coef)
                * (bi / b - 2.0 / a * part2)
                * np.log((zz + 2.414 * B_coef) / (zz - 0.414 * B_coef))
            )
            fugacity = np.exp(part1 + part3)

        elif rule == 2:
            _, gama = activity_fn(T, x, mixture.components, BIP_data)
            part1 = bi / b * (zz - 1.0) - np.log(zz - B_coef)
            sqrt2p1 = 1.0 + np.sqrt(2.0)
            sqrt2m1 = 1.0 - np.sqrt(2.0)
            part3 = (
                -1.0 / TwoSqrt2
                * (ai / (bi * R * T) - np.log(gama) / s1)
                * np.log((zz + sqrt2p1 * B_coef) / (zz + sqrt2m1 * B_coef))
            )
            fugacity = np.exp(part1 + part3)

        elif rule == 3:
            _, gama = activity_fn(T, x, mixture.components, BIP_data)
            q1 = -0.53
            sqrt2p1 = 1.0 + np.sqrt(2.0)
            sqrt2m1 = 1.0 - np.sqrt(2.0)
            logfi = (
                bi / b * (zz - 1.0)
                - np.log(zz - B_coef)
                - 1.0 / TwoSqrt2
                * (
                    ai / (bi * R * T)
                    + np.log(gama) / q1
                    + np.log(b / bi) / q1
                    + (bi / b - 1.0) / q1
                )
                * np.log((zz + sqrt2p1 * B_coef) / (zz + sqrt2m1 * B_coef))
            )
            fugacity = np.exp(logfi)

        elif rule == 4:
            _, gama = activity_fn(T, x, mixture.components, BIP_data)
            q1 = -0.4347
            q2 = -0.003654
            alphai = ai / (bi * R * T)
            alpha = a / (b * R * T)
            sqrt2p1 = 1.0 + np.sqrt(2.0)
            sqrt2m1 = 1.0 - np.sqrt(2.0)
            logfi = (
                bi / b * (zz - 1.0)
                - np.log(zz - B_coef)
                - 1.0 / TwoSqrt2
                * (
                    q1 * alphai
                    + q2 * (alpha ** 2 + alphai ** 2)
                    + np.log(gama)
                    + np.log(b / bi)
                    + bi / b
                    - 1.0
                )
                / (q1 + 2.0 * q2 * alpha)
                * np.log((zz + sqrt2p1 * B_coef) / (zz + sqrt2m1 * B_coef))
            )
            fugacity = np.exp(logfi)

    # Residual enthalpy (kij=0 approximation for dadT)
    sqrt_ai = np.sqrt(ai)
    sqrt_aci = np.sqrt(aci)
    dadT = (
        -float(np.dot(x, sqrt_ai))
        * np.dot(x, sqrt_aci * mi / np.sqrt(critical_temp))
        / np.sqrt(T)
    )

    TwoSqrt2 = 2.0 * np.sqrt(2.0)
    sqrt2p1 = 1.0 + np.sqrt(2.0)
    sqrt2m1 = 1.0 - np.sqrt(2.0)
    ln_term = np.log((zz + sqrt2p1 * B_coef) / (zz + sqrt2m1 * B_coef))
    HR = R * T * (zz - 1.0) + (T * dadT - a) / (b * TwoSqrt2) * ln_term

    props = None
    if compute_props:
        SR = R * np.log(zz - B_coef) + dadT / (b * TwoSqrt2) * ln_term
        GR = HR - T * SR
        VR = R * T * (zz - 1.0) / p

        dsqrtaidT = -sqrt_aci * mi / (2.0 * np.sqrt(critical_temp * T))
        d2sqrtaidT2 = sqrt_aci * mi / (4.0 * np.sqrt(critical_temp) * T ** 1.5)
        d2adT2 = (
            2.0 * float(np.dot(x, dsqrtaidT)) ** 2
            + 2.0 * float(np.dot(x, sqrt_ai)) * float(np.dot(x, d2sqrtaidT2))
        )

        Cv_R = -T * d2adT2 / (b * TwoSqrt2) * ln_term

        V_mol = zz * R * T / p
        den_pr = V_mol ** 2 + 2.0 * b * V_mol - b ** 2
        dPdT_V = R / (V_mol - b) - dadT / den_pr
        dPdV_T = -R * T / (V_mol - b) ** 2 + 2.0 * a * (V_mol + b) / den_pr ** 2
        Cp_R = Cv_R - T * dPdT_V ** 2 / dPdV_T - R

        props = ResidualProps(
            HR=float(HR), SR=float(SR), GR=float(GR), VR=float(VR),
            Cp_R=float(Cp_R), Cv_R=float(Cv_R),
        )

    return liquid_z, vapor_z, fugacity, float(HR), props

def SRKEOS(
    mixture: Mixture,
    thermo: ThermoModel,
    compute_props: bool = False,
) -> tuple[float, float, np.ndarray, float, ResidualProps | None]:
    """Soave-Redlich-Kwong (1972) Equation of State.

    Args:
        mixture: Mixture object with temperature, pressure, etc.
        thermo: ThermoModel with mixingrule, phase, fugacity_switch.
        compute_props: If True, compute and return full residual properties.

    Returns:
        Tuple of (liquid_z, vapor_z, fugacity, HR, props).
    """
    rule = thermo.mixingrule
    activity_fn = thermo.activity_model
    phase1 = thermo.phase
    fug_need = thermo.fugacity_switch

    critical_pres = np.array([c.Pc for c in mixture.components])
    critical_temp = np.array([c.Tc for c in mixture.components])
    omega = np.array([c.acentric_factor for c in mixture.components])
    BIP_data = mixture.bip
    x = mixture.mole_fraction
    p = mixture.pressure
    T = mixture.temperature
    N = len(critical_temp)

    fugacity = np.zeros(N)

    s1 = np.log(2.0)  # HV constant for SRK

    bi = 0.08664 * R * critical_temp / critical_pres
    aci = 0.42748 * (R * critical_temp) ** 2 / critical_pres
    mi = 0.48 + (1.574 - 0.176 * omega) * omega
    Tr = T / critical_temp
    alfai = 1.0 + mi * (1.0 - np.sqrt(Tr))
    alfa = alfai ** 2
    ai = aci * alfa

    Q = np.array([0.0, 0.0])
    if rule == 3:
        Q = np.array([-0.593, 0.0])
    elif rule == 4:
        Q = np.array([-0.478, -0.0047])

    a, b = mixing_rule(mixture, thermo, ai, bi, s1, Q)

    A_coef = a * p / (R * T) ** 2
    B_coef = b * p / (R * T)

    # SRK cubic: Z^3 - Z^2 + (A-B-B^2)Z - AB = 0
    poly_coef = np.array([1.0, -1.0, A_coef - B_coef - B_coef ** 2, -A_coef * B_coef])

    z_roots = np.roots(poly_coef)
    liquid_z, vapor_z = select_z_roots(z_roots, B_coef)

    if phase1 == 1:
        zz = liquid_z
    else:
        zz = vapor_z

    if fug_need == 1:
        if rule == 1:
            part1 = bi / b * (zz - 1.0) - np.log(zz - B_coef)
            a_matrix = np.sqrt(np.outer(ai, ai)) * (
                1.0 - BIP_data.EOScons - BIP_data.EOStdep * T
            )
            part2 = x @ a_matrix
            part3 = (
                A_coef / B_coef
                * (bi / b - 2.0 / a * part2)
                * np.log((zz + B_coef) / zz)
            )
            fugacity = np.exp(part1 + part3)

        elif rule == 2:
            _, gama = activity_fn(T, x, mixture.components, BIP_data)
            part1 = bi / b * (zz - 1.0) - np.log(zz - B_coef)
            part3 = -(ai / (bi * R * T) - np.log(gama) / s1) * np.log((zz + B_coef) / zz)
            fugacity = np.exp(part1 + part3)

        elif rule == 3:
            _, gama = activity_fn(T, x, mixture.components, BIP_data)
            q1 = -0.593
            logfi = (
                bi / b * (zz - 1.0)
                - np.log(zz - B_coef)
                - (
                    ai / (bi * R * T)
                    + np.log(gama) / q1
                    + np.log(b / bi) / q1
                    + (bi / b - 1.0) / q1
                )
                * np.log((zz + B_coef) / zz)
            )
            fugacity = np.exp(logfi)

        elif rule == 4:
            _, gama = activity_fn(T, x, mixture.components, BIP_data)
            q1 = -0.478
            q2 = -0.0047
            alphai = ai / (bi * R * T)
            alpha = a / (b * R * T)
            logfi = (
                bi / b * (zz - 1.0)
                - np.log(zz - B_coef)
                - (
                    q1 * alphai
                    + q2 * (alpha ** 2 + alphai ** 2)
                    + np.log(gama)
                    + np.log(b / bi)
                    + bi / b
                    - 1.0
                )
                / (q1 + 2.0 * q2 * alpha)
                * np.log((zz + B_coef) / zz)
            )
            fugacity = np.exp(logfi)

    # Residual enthalpy (kij=0 approximation)
    sqrt_ai = np.sqrt(ai)
    sqrt_aci = np.sqrt(aci)
    dadT = (
        -float(np.dot(x, sqrt_ai))
        * np.dot(x, sqrt_aci * mi / np.sqrt(critical_temp))
        / np.sqrt(T)
    )

    ln_term = np.log((zz + B_coef) / zz)
    HR = R * T * (zz - 1.0) + (T * dadT - a) / b * ln_term

    props = None
    if compute_props:
        SR = R * np.log(zz - B_coef) + dadT / b * ln_term
        GR = HR - T * SR
        VR = R * T * (zz - 1.0) / p

        dsqrtaidT = -sqrt_aci * mi / (2.0 * np.sqrt(critical_temp * T))
        d2sqrtaidT2 = sqrt_aci * mi / (4.0 * np.sqrt(critical_temp) * T ** 1.5)
        d2adT2 = (
            2.0 * float(np.dot(x, dsqrtaidT)) ** 2
            + 2.0 * float(np.dot(x, sqrt_ai)) * float(np.dot(x, d2sqrtaidT2))
        )

        Cv_R = -T * d2adT2 / b * ln_term

        V_mol = zz * R * T / p
        den_pr = V_mol ** 2 + 2.0 * b * V_mol - b ** 2
        dPdT_V = R / (V_mol - b) - dadT / den_pr
        dPdV_T = -R * T / (V_mol - b) ** 2 + 2.0 * a * (V_mol + b) / den_pr ** 2
        Cp_R = Cv_R - T * dPdT_V ** 2 / dPdV_T - R

        props = ResidualProps(
            HR=float(HR), SR=float(SR), GR=float(GR), VR=float(VR),
            Cp_R=float(Cp_R), Cv_R=float(Cv_R),
        )

    return liquid_z, vapor_z, fugacity, float(HR), props

def PR78EOS(
    mixture: Mixture,
    thermo: ThermoModel,
    compute_props: bool = False,
) -> tuple[float, float, np.ndarray, float, ResidualProps | None]:
    """Peng-Robinson EOS with 1978 alpha-function correction.

    For omega <= 0.491, uses the standard PR alpha function.
    For omega > 0.491, uses the 1978 modification (better for heavy compounds).

    Args:
        mixture: Mixture object.
        thermo: ThermoModel with mixingrule, phase, fugacity_switch.
        compute_props: If True, compute and return full residual properties.

    Returns:
        Tuple of (liquid_z, vapor_z, fugacity, HR, props).
    """
    rule = thermo.mixingrule
    activity_fn = thermo.activity_model
    phase1 = thermo.phase
    fug_need = thermo.fugacity_switch

    critical_pres = np.array([c.Pc for c in mixture.components])
    critical_temp = np.array([c.Tc for c in mixture.components])
    omega = np.array([c.acentric_factor for c in mixture.components])
    BIP_data = mixture.bip
    x = mixture.mole_fraction
    p = mixture.pressure
    T = mixture.temperature
    N = len(critical_temp)

    fugacity = np.zeros(N)

    s1 = 0.623225

    bi = 0.077796 * R * critical_temp / critical_pres
    aci = 0.457235 * (R * critical_temp) ** 2 / critical_pres

    mi = 0.37646 + (1.54226 - 0.26992 * omega) * omega
    high_omega = omega > 0.491
    mi[high_omega] = (
        0.379642
        + 1.48503 * omega[high_omega]
        - 0.164423 * omega[high_omega] ** 2
        + 0.016666 * omega[high_omega] ** 3
    )

    Tr = T / critical_temp
    alfai = 1.0 + mi * (1.0 - np.sqrt(Tr))
    alfa = alfai ** 2
    ai = aci * alfa

    Q = np.array([0.0, 0.0])
    if rule == 3:
        Q = np.array([-0.53, 0.0])
    elif rule == 4:
        Q = np.array([-0.4347, -0.003654])

    a, b = mixing_rule(mixture, thermo, ai, bi, s1, Q)

    A_coef = a * p / (R * T) ** 2
    B_coef = b * p / (R * T)

    poly_coef = np.array([
        1.0,
        -1.0 + B_coef,
        A_coef - B_coef * (2.0 + 3.0 * B_coef),
        -B_coef * (A_coef - B_coef * (1.0 + B_coef)),
    ])

    z_roots = np.roots(poly_coef)
    liquid_z, vapor_z = select_z_roots(z_roots, B_coef)

    if phase1 == 1:
        zz = liquid_z
    else:
        zz = vapor_z

    if fug_need == 1:
        TwoSqrt2 = 2.0 * np.sqrt(2.0)
        sqrt2p1 = 1.0 + np.sqrt(2.0)
        sqrt2m1 = 1.0 - np.sqrt(2.0)

        if rule == 1:
            part1 = bi / b * (zz - 1.0) - np.log(zz - B_coef)
            a_matrix = np.sqrt(np.outer(ai, ai)) * (
                1.0 - BIP_data.EOScons - BIP_data.EOStdep * T
            )
            part2 = x @ a_matrix
            part3 = (
                A_coef / (2.828 * B_coef)
                * (bi / b - 2.0 / a * part2)
                * np.log((zz + 2.414 * B_coef) / (zz - 0.414 * B_coef))
            )
            fugacity = np.exp(part1 + part3)

        elif rule == 2:
            _, gama = activity_fn(T, x, mixture.components, BIP_data)
            part1 = bi / b * (zz - 1.0) - np.log(zz - B_coef)
            part3 = (
                -1.0 / TwoSqrt2
                * (ai / (bi * R * T) - np.log(gama) / s1)
                * np.log((zz + sqrt2p1 * B_coef) / (zz + sqrt2m1 * B_coef))
            )
            fugacity = np.exp(part1 + part3)

        elif rule == 3:
            _, gama = activity_fn(T, x, mixture.components, BIP_data)
            q1 = -0.53
            logfi = (
                bi / b * (zz - 1.0)
                - np.log(zz - B_coef)
                - 1.0 / TwoSqrt2
                * (
                    ai / (bi * R * T)
                    + np.log(gama) / q1
                    + np.log(b / bi) / q1
                    + (bi / b - 1.0) / q1
                )
                * np.log((zz + sqrt2p1 * B_coef) / (zz + sqrt2m1 * B_coef))
            )
            fugacity = np.exp(logfi)

        elif rule == 4:
            _, gama = activity_fn(T, x, mixture.components, BIP_data)
            q1 = -0.4347
            q2 = -0.003654
            alphai = ai / (bi * R * T)
            alpha = a / (b * R * T)
            logfi = (
                bi / b * (zz - 1.0)
                - np.log(zz - B_coef)
                - 1.0 / TwoSqrt2
                * (
                    q1 * alphai
                    + q2 * (alpha ** 2 + alphai ** 2)
                    + np.log(gama)
                    + np.log(b / bi)
                    + bi / b
                    - 1.0
                )
                / (q1 + 2.0 * q2 * alpha)
                * np.log((zz + sqrt2p1 * B_coef) / (zz + sqrt2m1 * B_coef))
            )
            fugacity = np.exp(logfi)

    sqrt_ai = np.sqrt(ai)
    sqrt_aci = np.sqrt(aci)
    dadT = (
        -float(np.dot(x, sqrt_ai))
        * np.dot(x, sqrt_aci * mi / np.sqrt(critical_temp))
        / np.sqrt(T)
    )

    TwoSqrt2 = 2.0 * np.sqrt(2.0)
    sqrt2p1 = 1.0 + np.sqrt(2.0)
    sqrt2m1 = 1.0 - np.sqrt(2.0)
    ln_term = np.log((zz + sqrt2p1 * B_coef) / (zz + sqrt2m1 * B_coef))
    HR = R * T * (zz - 1.0) + (T * dadT - a) / (b * TwoSqrt2) * ln_term

    props = None
    if compute_props:
        SR = R * np.log(zz - B_coef) + dadT / (b * TwoSqrt2) * ln_term
        GR = HR - T * SR
        VR = R * T * (zz - 1.0) / p

        dsqrtaidT = -sqrt_aci * mi / (2.0 * np.sqrt(critical_temp * T))
        d2sqrtaidT2 = sqrt_aci * mi / (4.0 * np.sqrt(critical_temp) * T ** 1.5)
        d2adT2 = (
            2.0 * float(np.dot(x, dsqrtaidT)) ** 2
            + 2.0 * float(np.dot(x, sqrt_ai)) * float(np.dot(x, d2sqrtaidT2))
        )

        Cv_R = -T * d2adT2 / (b * TwoSqrt2) * ln_term

        V_mol = zz * R * T / p
        den_pr = V_mol ** 2 + 2.0 * b * V_mol - b ** 2
        dPdT_V = R / (V_mol - b) - dadT / den_pr
        dPdV_T = -R * T / (V_mol - b) ** 2 + 2.0 * a * (V_mol + b) / den_pr ** 2
        Cp_R = Cv_R - T * dPdT_V ** 2 / dPdV_T - R

        props = ResidualProps(
            HR=float(HR), SR=float(SR), GR=float(GR), VR=float(VR),
            Cp_R=float(Cp_R), Cv_R=float(Cv_R),
        )

    return liquid_z, vapor_z, fugacity, float(HR), props
