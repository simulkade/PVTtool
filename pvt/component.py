"""Component — pure-component thermodynamic properties.

Each Component holds critical properties, correlations for vapor pressure,
enthalpy of vaporization, liquid and ideal-gas heat capacity, and formation
data. Components are immutable value objects.
"""

from dataclasses import dataclass, field

import numpy as np

from ._data import get_component_data, get_components_data
from ._equations import EQUATION_FUNCTIONS


@dataclass(frozen=True)
class Component:
    """Pure-component thermodynamic properties container.

    All properties are in SI units: temperature [K], pressure [Pa], energy [J],
    heat capacity [J/(mol·K)], volume [m³/mol].

    Attributes:
        name: Component name (e.g. 'Methane').
        formula: Chemical formula (e.g. 'CH4').
        MW: Molecular weight [kg/mol].
        Tc: Critical temperature [K].
        Pc: Critical pressure [Pa].
        Vc: Critical volume [m³/mol].
        Zc: Critical compressibility factor.
        acentric_factor: Pitzer acentric factor.
        Psat_eq: Vapor pressure equation type identifier.
        Psat_coefs: Vapor pressure coefficients (1D array).
        Psat_Trange: Valid temperature range [Tmin, Tmax] for Psat.
        Psat_range: Pressure range corresponding to Psat_Trange.
        dh_vap_eq: Enthalpy of vaporization equation type.
        dh_vap_coefs: dh_vap coefficients (1D array).
        dh_vap_Trange: Valid T range [Tmin, Tmax] for dh_vap.
        dh_vap_range: dh_vap values at Trange bounds.
        cp_liq_eq: Liquid heat capacity equation type.
        cp_liq_coefs: cp_liq coefficients (1D array).
        cp_liq_Trange: Valid T range for cp_liq.
        cp_liq_range: cp_liq values at Trange bounds.
        cp_ig_eq: Ideal-gas heat capacity equation type.
        cp_ig_coefs: cp_ig coefficients (1D array).
        cp_ig_Trange: Valid T range for cp_ig.
        cp_ig_range: cp_ig values at Trange bounds.
        dhf_ig: Ideal-gas enthalpy of formation [J/mol].
        dgf_ig: Ideal-gas Gibbs free energy of formation [J/mol].
        ds_ig: Ideal-gas absolute entropy [J/(mol·K)].
        dh_comb: Enthalpy of combustion [J/mol].
    """

    name: str
    formula: str
    MW: float
    Tc: float
    Pc: float
    Vc: float
    Zc: float
    acentric_factor: float
    Psat_eq: str
    Psat_coefs: np.ndarray = field(default_factory=lambda: np.zeros(5))
    Psat_Trange: np.ndarray = field(default_factory=lambda: np.zeros(2))
    Psat_range: np.ndarray = field(default_factory=lambda: np.zeros(2))
    dh_vap_eq: str = ""
    dh_vap_coefs: np.ndarray = field(default_factory=lambda: np.zeros(4))
    dh_vap_Trange: np.ndarray = field(default_factory=lambda: np.zeros(2))
    dh_vap_range: np.ndarray = field(default_factory=lambda: np.zeros(2))
    cp_liq_eq: str = ""
    cp_liq_coefs: np.ndarray = field(default_factory=lambda: np.zeros(5))
    cp_liq_Trange: np.ndarray = field(default_factory=lambda: np.zeros(2))
    cp_liq_range: np.ndarray = field(default_factory=lambda: np.zeros(2))
    cp_ig_eq: str = ""
    cp_ig_coefs: np.ndarray = field(default_factory=lambda: np.zeros(5))
    cp_ig_Trange: np.ndarray = field(default_factory=lambda: np.zeros(2))
    cp_ig_range: np.ndarray = field(default_factory=lambda: np.zeros(2))
    dhf_ig: float = 0.0
    dgf_ig: float = 0.0
    ds_ig: float = 0.0
    dh_comb: float = 0.0

    def vapor_pressure(self, T: float) -> float:
        """Compute vapor pressure at temperature T [K].

        Returns p_sat in [Pa].
        """
        fn = EQUATION_FUNCTIONS[self.Psat_eq]
        return fn(T, self.Psat_coefs)

    def dh_vap(self, T: float) -> float:
        """Compute enthalpy of vaporization at temperature T [K].

        Returns dh_vap in [J/mol].
        """
        fn = EQUATION_FUNCTIONS[self.dh_vap_eq]
        return fn(T / self.Tc, self.dh_vap_coefs)

    def cp_liq(self, T: float) -> float:
        """Compute liquid heat capacity at temperature T [K].

        Returns cp in [J/(mol·K)].
        """
        fn = EQUATION_FUNCTIONS[self.cp_liq_eq]
        if "dippr100" in self.cp_liq_eq:
            return fn(T / self.Tc, self.cp_liq_coefs)
        return fn(T, self.cp_liq_coefs)

    def cp_ig(self, T: float) -> float:
        """Compute ideal-gas heat capacity at temperature T [K].

        Returns cp in [J/(mol·K)].
        """
        fn = EQUATION_FUNCTIONS[self.cp_ig_eq]
        return fn(T, self.cp_ig_coefs)

    @staticmethod
    def from_database(name_or_formula: str) -> "Component":
        """Load a single component from the built-in database.

        Args:
            name_or_formula: Component name (e.g. 'Methane') or formula (e.g. 'CH4').

        Returns:
            Component instance.

        Raises:
            KeyError: If the component is not found in the database.
        """
        data = get_component_data(name_or_formula)
        if data is None:
            raise KeyError(f"Component '{name_or_formula}' not found in database.")
        return Component._from_dict(data)

    @staticmethod
    def from_database_array(names_or_formulas: list[str]) -> tuple[list["Component"], list[int]]:
        """Load multiple components from the built-in database.

        Args:
            names_or_formulas: List of component names or formulas.

        Returns:
            Tuple of (components, not_found_indices).
        """
        data_list, not_found = get_components_data(names_or_formulas)
        components = [Component._from_dict(d) for d in data_list]
        return components, not_found

    @staticmethod
    def _from_dict(data: dict) -> "Component":
        """Create a Component from a database dictionary."""
        return Component(
            name=data["name"],
            formula=data["formula"],
            MW=float(data["MW"]),
            Tc=float(data["Tc"]),
            Pc=float(data["Pc"]),
            Vc=float(data["Vc"]),
            Zc=float(data["Zc"]),
            acentric_factor=float(data["acentric_factor"]),
            Psat_eq=data["Psat_eq"],
            Psat_coefs=np.array(data["Psat_coefs"]),
            Psat_Trange=np.array(data["Psat_Trange"]),
            Psat_range=np.array(data["Psat_range"]),
            dh_vap_eq=data["dh_vap_eq"],
            dh_vap_coefs=np.array(data["dh_vap_coefs"]),
            dh_vap_Trange=np.array(data["dh_vap_Trange"]),
            dh_vap_range=np.array(data["dh_vap_range"]),
            cp_liq_eq=data["cp_liq_eq"],
            cp_liq_coefs=np.array(data["cp_liq_coefs"]),
            cp_liq_Trange=np.array(data["cp_liq_Trange"]),
            cp_liq_range=np.array(data["cp_liq_range"]),
            cp_ig_eq=data["cp_ig_eq"],
            cp_ig_coefs=np.array(data["cp_ig_coefs"]),
            cp_ig_Trange=np.array(data["cp_ig_Trange"]),
            cp_ig_range=np.array(data["cp_ig_range"]),
            dhf_ig=float(data["dhf_ig"]),
            dgf_ig=float(data["dgf_ig"]),
            ds_ig=float(data["ds_ig"]),
            dh_comb=float(data["dh_comb"]),
        )
