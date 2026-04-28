"""One-shot script: convert PureData/puredata.mat to src/pvttool/data/puredata.json."""

import json
import pathlib
import sys

import numpy as np
import scipy.io


def _to_python(obj):
    """Recursively convert numpy scalars/arrays to plain Python types."""
    if isinstance(obj, np.ndarray):
        if obj.ndim == 0:
            return _to_python(obj.item())
        return [_to_python(v) for v in obj.tolist()]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.complexfloating,)):
        return float(obj.real)
    return obj


def main():
    repo = pathlib.Path(__file__).parent.parent
    mat_path = repo / "PureData" / "puredata.mat"
    out_path = repo / "src" / "pvttool" / "data" / "puredata.json"

    if not mat_path.exists():
        print(f"ERROR: {mat_path} not found", file=sys.stderr)
        sys.exit(1)

    raw = scipy.io.loadmat(str(mat_path), squeeze_me=True)

    def arr(key):
        v = raw[key]
        if isinstance(v, np.ndarray):
            return v.tolist()
        return v

    def scalar_arr(key):
        v = raw[key]
        if isinstance(v, np.ndarray):
            return [float(x) for x in v.flat]
        return [float(v)]

    def str_arr(key):
        v = raw[key]
        if isinstance(v, np.ndarray) and v.dtype.kind in ("U", "S", "O"):
            return [str(s).strip() for s in v.flat]
        return [str(v).strip()]

    def matrix(key):
        v = raw[key]
        if isinstance(v, np.ndarray):
            if v.ndim == 1:
                return [[float(x)] for x in v]
            return [[float(v[i, j]) for j in range(v.shape[1])] for i in range(v.shape[0])]
        return [[float(v)]]

    def scalar_int_arr(key):
        v = raw[key]
        if isinstance(v, np.ndarray):
            return [int(x) for x in v.flat]
        return [int(v)]

    db = {
        "component_name": str_arr("component_name"),
        "component_formula": str_arr("component_formula"),
        "MW": scalar_arr("Molecular_weight_data"),
        "Tc": scalar_arr("critical_temperature_data"),
        "Pc": scalar_arr("critical_pressure_data"),
        "Vc": scalar_arr("critical_volume_data"),
        "Zc": scalar_arr("critical_compressibility_factor"),
        "acentric_factor": scalar_arr("acentric_factor"),
        # vapor pressure
        "vapor_pressure_coefs": matrix("vapor_pressure_coefs"),
        "vapor_pressure_T_range": matrix("vapor_pressure_T_range"),
        "vapor_pressure_p_range": matrix("vapor_pressure_p_range"),
        # dh vaporization
        "dh_vaporization_coefs": matrix("dh_vaporization_coefs"),
        "dh_vaporization_T_range": matrix("dh_vaporization_T_range"),
        "dh_vaporization_dhv_range": matrix("dh_vaporization_dhv_range"),
        # liquid heat capacity
        "liquid_heat_capacity_equation_num": scalar_int_arr("liquid_heat_capacity_equation_num"),
        "liquid_heat_capacity_coef": matrix("liquid_heat_capacity_coef"),
        "liquid_heat_capacity_T_range": matrix("liquid_heat_capacity_T_range"),
        "liquid_heat_capacity_cp_range": matrix("liquid_heat_capacity_cp_range"),
        # ideal gas heat capacity
        "ideal_gas_heat_capacity_equation_num": scalar_int_arr("ideal_gas_heat_capacity_equation_num"),
        "ideal_gas_heat_capacity_coefs": matrix("ideal_gas_heat_capacity_coefs"),
        "ideal_gas_heat_capacity_T_range": matrix("ideal_gas_heat_capacity_T_range"),
        "ideal_gas_heat_capacity_cp_range": matrix("ideal_gas_heat_capacity_cp_range"),
        # formation and combustion
        "dh_formation_ig": scalar_arr("dh_formation_ig"),
        "dg_formation_ig": scalar_arr("dg_formation_ig"),
        "ds_ideal_gas": scalar_arr("ds_ideal_gas"),
        "dh_combustion": scalar_arr("dh_combustion"),
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(db, f, indent=2)

    n = len(db["component_name"])
    print(f"Wrote {n} components to {out_path}")


if __name__ == "__main__":
    main()
