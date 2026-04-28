"""Pure-component thermodynamic properties."""

from __future__ import annotations

import functools
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# Correlation functions (translated directly from MATLAB function handles)
# ---------------------------------------------------------------------------

def _vapor_pressure(T: float, coef: list[float]) -> float:
    """Extended Antoine: Psat = exp(c0 + c1/T + c2*ln(T) + c3*T^c4) [Pa]."""
    return math.exp(coef[0] + coef[1] / T + coef[2] * math.log(T) + coef[3] * T ** coef[4])


def _dh_vap(Tr: float, coef: list[float]) -> float:
    """Watson-type: dh = c0*(1-Tr)^(c1+c2*Tr+c3*Tr^2) / 1000 [J/mol]."""
    exp = coef[1] + coef[2] * Tr + coef[3] * Tr ** 2
    return coef[0] * (1.0 - Tr) ** exp / 1000.0


def _cp_liq_1(Tr: float, coef: list[float]) -> float:
    """Polynomial in Tr: (c0 + c1*Tr + c2*Tr^2 + c3*Tr^3 + c4*Tr^4) / 1000 [J/(mol*K)]."""
    return (coef[0] + coef[1] * Tr + coef[2] * Tr**2 + coef[3] * Tr**3 + coef[4] * Tr**4) / 1000.0


def _cp_liq_2(Tr: float, coef: list[float]) -> float:
    """Riedel form: (c0^2/(1-Tr) + c1 - 2*c0*c2*(1-Tr) - ...) / 1000 [J/(mol*K)]."""
    dt = 1.0 - Tr
    val = (
        coef[0] ** 2 / dt
        + coef[1]
        - 2.0 * coef[0] * coef[2] * dt
        - coef[0] * coef[3] * dt**2
        - coef[2] ** (2.0 / 3.0) * dt**3
        - coef[2] * coef[3] / 2.0 * dt**4
        - coef[3] ** (2.0 / 5.0) * dt**5
    )
    return val / 1000.0


def _cp_ig_1(T: float, coef: list[float]) -> float:
    """DIPPR 107 hyperbolic: (c0 + c1*(c2/T/sinh(c2/T))^2 + c3*(c4/T/cosh(c4/T))^2) / 1000."""
    t1 = coef[2] / T
    t2 = coef[4] / T
    sinh_t1 = math.sinh(t1)
    cosh_t2 = math.cosh(t2)
    val = coef[0] + coef[1] * (t1 / sinh_t1) ** 2 + coef[3] * (t2 / cosh_t2) ** 2
    return val / 1000.0


def _cp_ig_2(T: float, coef: list[float]) -> float:
    """Polynomial in T: (c0 + c1*T + c2*T^2 + c3*T^3 + c4*T^4) / 1000 [J/(mol*K)]."""
    return (coef[0] + coef[1] * T + coef[2] * T**2 + coef[3] * T**3 + coef[4] * T**4) / 1000.0


def _cp_ig_3(T: float, coef: list[float]) -> float:
    """Log form: (c0 + c1*ln(T) + c2/T + c3*T) / 1000 [J/(mol*K)]."""
    return (coef[0] + coef[1] * math.log(T) + coef[2] / T + coef[3] * T) / 1000.0


# ---------------------------------------------------------------------------
# Database loader
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=1)
def _load_database() -> dict:
    """Load puredata.json (cached after first call)."""
    db_path = Path(__file__).parent.parent / "data" / "puredata.json"
    with open(db_path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Component dataclass
# ---------------------------------------------------------------------------

@dataclass
class Component:
    """Pure-component thermodynamic properties.

    Parameters
    ----------
    name : str
        Full component name (e.g. "Methane").
    formula : str
        Chemical formula (e.g. "CH4").
    MW : float
        Molar mass [kg/mol].
    Tc : float
        Critical temperature [K].
    Pc : float
        Critical pressure [Pa].
    Vc : float
        Critical molar volume [m³/mol].
    Zc : float
        Critical compressibility factor [-].
    acentric_factor : float
        Pitzer acentric factor [-].
    psat_coefs : list[float]
        Coefficients for the extended-Antoine vapor-pressure equation (5 values).
    psat_T_range : list[float]
        Valid temperature range [T_min, T_max] [K] for vapor pressure.
    dh_vap_coefs : list[float]
        Watson-type dHvap correlation coefficients (4 values).
    dh_vap_T_range : list[float]
        Valid Tr range for dHvap.
    cp_liq_eq_num : int
        Which liquid Cp equation to use (1 or 2).
    cp_liq_coefs : list[float]
        Liquid heat-capacity correlation coefficients (5 values).
    cp_liq_T_range : list[float]
        Valid Tr range for liquid Cp.
    cp_ig_eq_num : int
        Which ideal-gas Cp equation to use (1, 2, or 3).
    cp_ig_coefs : list[float]
        Ideal-gas heat-capacity correlation coefficients (5 values).
    cp_ig_T_range : list[float]
        Valid T range [K] for ideal-gas Cp.
    dhf_ig : float
        Standard ideal-gas enthalpy of formation [J/mol].
    dgf_ig : float
        Standard ideal-gas Gibbs energy of formation [J/mol].
    ds_ig : float
        Standard ideal-gas entropy [J/(mol*K)].
    dh_comb : float
        Standard enthalpy of combustion [J/mol].
    """

    name: str
    formula: str
    MW: float
    Tc: float
    Pc: float
    Vc: float
    Zc: float
    acentric_factor: float
    psat_coefs: list
    psat_T_range: list
    dh_vap_coefs: list
    dh_vap_T_range: list
    cp_liq_eq_num: int
    cp_liq_coefs: list
    cp_liq_T_range: list
    cp_ig_eq_num: int
    cp_ig_coefs: list
    cp_ig_T_range: list
    dhf_ig: float
    dgf_ig: float
    ds_ig: float
    dh_comb: float

    # Extra fields for activity models (UNIQUAC)
    uniquacR: float = 0.0
    uniquacQ: float = 0.0

    def vapor_pressure(self, T: float) -> float:
        """Return saturation pressure [Pa] at temperature T [K].

        Parameters
        ----------
        T : float
            Temperature [K].

        Returns
        -------
        float
            Saturation pressure [Pa].
        """
        return _vapor_pressure(T, self.psat_coefs)

    def dh_vap(self, T: float) -> float:
        """Return molar enthalpy of vaporization [J/mol] at T [K].

        Parameters
        ----------
        T : float
            Temperature [K].

        Returns
        -------
        float
            Enthalpy of vaporization [J/mol].
        """
        Tr = T / self.Tc
        return _dh_vap(Tr, self.dh_vap_coefs)

    def cp_liq(self, T: float) -> float:
        """Return liquid molar heat capacity [J/(mol*K)] at T [K].

        Parameters
        ----------
        T : float
            Temperature [K].

        Returns
        -------
        float
            Liquid heat capacity [J/(mol*K)].
        """
        Tr = T / self.Tc
        if self.cp_liq_eq_num == 1:
            return _cp_liq_1(Tr, self.cp_liq_coefs)
        return _cp_liq_2(Tr, self.cp_liq_coefs)

    def cp_ig(self, T: float) -> float:
        """Return ideal-gas molar heat capacity [J/(mol*K)] at T [K].

        Parameters
        ----------
        T : float
            Temperature [K].

        Returns
        -------
        float
            Ideal-gas heat capacity [J/(mol*K)].
        """
        if self.cp_ig_eq_num == 1:
            return _cp_ig_1(T, self.cp_ig_coefs)
        elif self.cp_ig_eq_num == 2:
            return _cp_ig_2(T, self.cp_ig_coefs)
        return _cp_ig_3(T, self.cp_ig_coefs)

    @classmethod
    def from_database(cls, name: str) -> "Component":
        """Load a single component from the built-in database.

        Parameters
        ----------
        name : str
            Component name or chemical formula (case-insensitive).

        Returns
        -------
        Component
            The matching component.

        Raises
        ------
        KeyError
            If the component is not found in the database.
        """
        comps = cls.from_database_array([name])
        if not comps:
            raise KeyError(f"Component '{name}' not found in database.")
        return comps[0]

    @classmethod
    def from_database_array(cls, names: list[str]) -> list["Component"]:
        """Load multiple components from the built-in database.

        Parameters
        ----------
        names : list[str]
            Component names or chemical formulae (case-insensitive).

        Returns
        -------
        list[Component]
            List of matched components (in the order of `names`).

        Raises
        ------
        KeyError
            If any name is not found.
        """
        db = _load_database()
        db_names = [s.strip().lower() for s in db["component_name"]]
        db_formulas = [s.strip().lower() for s in db["component_formula"]]
        components = []
        for name in names:
            key = name.strip().lower()
            idx = None
            if key in db_formulas:
                idx = db_formulas.index(key)
            elif key in db_names:
                idx = db_names.index(key)
            if idx is None:
                raise KeyError(f"Component '{name}' not found in database.")
            components.append(cls._from_db_index(db, idx))
        return components

    @classmethod
    def _from_db_index(cls, db: dict, i: int) -> "Component":
        """Construct a Component from database dict at row index i."""
        return cls(
            name=db["component_name"][i],
            formula=db["component_formula"][i],
            MW=db["MW"][i],
            Tc=db["Tc"][i],
            Pc=db["Pc"][i],
            Vc=db["Vc"][i],
            Zc=db["Zc"][i],
            acentric_factor=db["acentric_factor"][i],
            psat_coefs=db["vapor_pressure_coefs"][i],
            psat_T_range=db["vapor_pressure_T_range"][i],
            dh_vap_coefs=db["dh_vaporization_coefs"][i],
            dh_vap_T_range=db["dh_vaporization_T_range"][i],
            cp_liq_eq_num=db["liquid_heat_capacity_equation_num"][i],
            cp_liq_coefs=db["liquid_heat_capacity_coef"][i],
            cp_liq_T_range=db["liquid_heat_capacity_T_range"][i],
            cp_ig_eq_num=db["ideal_gas_heat_capacity_equation_num"][i],
            cp_ig_coefs=db["ideal_gas_heat_capacity_coefs"][i],
            cp_ig_T_range=db["ideal_gas_heat_capacity_T_range"][i],
            dhf_ig=db["dh_formation_ig"][i],
            dgf_ig=db["dg_formation_ig"][i],
            ds_ig=db["ds_ideal_gas"][i],
            dh_comb=db["dh_combustion"][i],
        )
