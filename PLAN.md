# PVTtool Python Port Plan

## Context

PVTtool is a mature MATLAB/Octave thermodynamics library for PVT (Pressure-Volume-Temperature)
calculations. It implements cubic equations of state (PR, SRK, PR78), phase equilibrium flash
(VLE/LLE), stability testing, and residual property computation for multicomponent mixtures.

This plan converts the 61-file MATLAB package to a modern, installable Python package managed
with `uv`. The user-facing API is kept as close as possible to the MATLAB version, translated
to Python naming conventions (snake_case). All classes and public functions receive NumPy-style
docstrings and a full pytest test suite is included.

---

## Branch

```bash
git checkout -b pvt-claudecode-sonnet
```

---

## Target Package Layout

```
pvttool/                         ← repo root (alongside existing MATLAB files)
├── pyproject.toml
├── PLAN.md                      ← this file
├── scripts/
│   └── convert_puredata.py      ← one-shot: puredata.mat → puredata.json
├── src/
│   └── pvttool/
│       ├── __init__.py          ← public API re-exports
│       ├── classes/
│       │   ├── __init__.py
│       │   ├── component.py     ← Component dataclass + from_database()
│       │   ├── bip.py           ← BIP dataclass (n×n NumPy matrices)
│       │   ├── mixture.py       ← Mixture dataclass
│       │   ├── thermo_model.py  ← ThermoModel config dataclass
│       │   └── flash_options.py ← FlashOptions config dataclass
│       ├── eos/
│       │   ├── __init__.py
│       │   ├── pr.py            ← preos()
│       │   ├── srk.py           ← srkeos()
│       │   ├── pr78.py          ← pr78eos()
│       │   └── _roots.py        ← select_z_roots() (private)
│       ├── flash/
│       │   ├── __init__.py
│       │   ├── vle.py           ← vle_flash(), vle_flash_negative()
│       │   ├── lle.py           ← lle_flash()
│       │   ├── saturation.py    ← bubble_pressure/temperature, dew_pressure/temperature
│       │   ├── stability.py     ← stability_test(), stability_lle_test()
│       │   ├── _kvalue.py       ← kvalue(), kvalue_lle(), kval_estimate()
│       │   ├── _rachford_rice.py← rachford_rice_nr(), mass_bal_func(), xy_calc()
│       │   └── _fugacity.py     ← fugacity()
│       ├── activity/
│       │   ├── __init__.py
│       │   ├── nrtl.py          ← nrtl()
│       │   ├── wilson.py        ← wilson()
│       │   └── uniquac.py       ← uniquac()
│       ├── mixing_rules.py      ← mixing_rule()
│       ├── tools.py             ← add_components(), add_mixture(), add_thermo(), zero_bip()
│       ├── auxiliary.py         ← wilson_correlation(), normalize()
│       └── data/
│           └── puredata.json    ← converted from PureData/puredata.mat
└── tests/
    ├── conftest.py              ← shared fixtures (methane, methanol/water mixture, etc.)
    ├── test_classes.py
    ├── test_eos.py
    ├── test_flash_vle.py
    ├── test_flash_lle.py
    ├── test_saturation.py
    ├── test_stability.py
    └── test_residual_props.py
```

---

## Step 1: pyproject.toml

```toml
[project]
name = "pvttool"
version = "0.1.0"
requires-python = ">=3.11"
dependencies = [
    "numpy>=1.26",
    "scipy>=1.12",
]

[project.optional-dependencies]
dev = ["pytest>=8.0", "pytest-cov"]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[tool.hatch.build.targets.wheel]
packages = ["src/pvttool"]
```

Install: `uv pip install -e ".[dev]"`

---

## Step 2: Data Migration

`PureData/puredata.mat` is a MATLAB binary. Convert it once with a helper script:

```bash
python scripts/convert_puredata.py
# reads PureData/puredata.mat via scipy.io.loadmat
# writes src/pvttool/data/puredata.json
```

`Component.from_database()` lazy-loads the JSON (cached with `functools.lru_cache`).

---

## Step 3: Class Design

All classes use `@dataclass` for value semantics, matching MATLAB's value class behaviour.

### `classes/component.py` — `Component`

```python
@dataclass
class Component:
    name: str
    formula: str
    MW: float            # kg/mol
    Tc: float            # K
    Pc: float            # Pa
    Vc: float            # m³/mol
    Zc: float
    acentric_factor: float
    # correlation equation type + coefficients + T range for:
    # vapor_pressure, dh_vap, cp_liq, cp_ig, dhf_ig, dgf_ig, ds_ig, dh_comb

    def vapor_pressure(self, T: float) -> float: ...
    def dh_vap(self, T: float) -> float: ...
    def cp_liq(self, T: float) -> float: ...
    def cp_ig(self, T: float) -> float: ...

    @classmethod
    def from_database(cls, name: str) -> "Component": ...
    @classmethod
    def from_database_array(cls, names: list[str]) -> list["Component"]: ...
```

### `classes/bip.py` — `BIP`

```python
@dataclass
class BIP:
    n: int
    # all n×n arrays, initialised to zeros in __post_init__
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
    uniquac_r: np.ndarray = field(init=False)   # shape (n,)
    uniquac_q: np.ndarray = field(init=False)   # shape (n,)
```

### `classes/mixture.py` — `Mixture`

```python
@dataclass
class Mixture:
    components: list[Component]
    temperature: float          # K
    pressure: float             # Pa
    mole_fraction: np.ndarray   # shape (n,), defaults to equimolar
    bip: BIP = field(default_factory=...)
```

### `classes/thermo_model.py` — `ThermoModel`

```python
@dataclass
class ThermoModel:
    eos: Callable = field(default_factory=lambda: preos)
    activity_model: Callable = field(default_factory=lambda: nrtl)
    mixing_rule: int = 1     # 1=vdW, 2=HV, 3=MHV1, 4=MHV2
    phase: int = 1           # 1=liquid Z-root, 2=vapor Z-root
    fugacity_switch: int = 1
```

### `classes/flash_options.py` — `FlashOptions`

```python
@dataclass
class FlashOptions:
    accuracy: float = 1e-7
    iteration: int = 100
    trivial_solution_max_error: float = 1e-5
    convergence_max_error: float = 1e-10
    max_iteration: int = 50
```

---

## Step 4: EOS Functions

**Signature (all three EOS):**
```python
def preos(mixture: Mixture, thermo: ThermoModel) \
    -> tuple[float, float, np.ndarray, float, dict]:
    ...
    return z_liquid, z_vapor, fugacity_coefficients, HR, props
```

`props` dict keys: `HR`, `SR`, `GR`, `VR`, `Cp_R`, `Cv_R`.

- Cubic Z-roots via `np.roots()`.
- Root selection via `eos/_roots.py::select_z_roots(z_roots, B)`:
  - Discard roots with `|imag(Z)| / max(|real(Z)|, 1) > 1e-6`.
  - Discard roots with `real(Z) <= B`.
  - Return min and max of remaining physical roots.
- `mixing_rule()` called before fugacity computation.
- All three EOS share the same `_roots.select_z_roots` helper.

---

## Step 5: Mixing Rules

`mixing_rules.py::mixing_rule(mixture, thermo, ai, bi, s1, Q) -> tuple[float, float]`

Supports rules 1–4 (van der Waals, Huron-Vidal, MHV1, MHV2). Rule 5
(Wong-Sandler) raises `NotImplementedError`.

---

## Step 6: Activity Models

**Common signature:**
```python
def nrtl(T: float, x: np.ndarray, components: list[Component], bip: BIP) \
    -> tuple[float, np.ndarray]:
    ...
    return g_ert, gamma   # dimensionless excess Gibbs, activity coefficients
```

Same pattern for `wilson()` and `uniquac()`.

---

## Step 7: Flash Functions

### VLE Flash (`flash/vle.py`)

```python
def vle_flash(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Returns (vapor_y, liquid_x, vapor_fraction)."""
```

Algorithm: Rachford-Rice solved via Newton-Raphson (`_rachford_rice.py`),
with GDEM acceleration every 5 iterations.

### LLE Flash (`flash/lle.py`)

```python
def lle_flash(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Returns (liquid2_y, liquid1_x, liquid2_fraction)."""
```

Both phases evaluated with `thermo.phase = 1` (liquid Z-root).

### Saturation (`flash/saturation.py`)

```python
def bubble_pressure(mixture, thermo, options=None) -> tuple[float, np.ndarray, bool]: ...
def bubble_temperature(mixture, thermo, options=None) -> tuple[float, np.ndarray, bool]: ...
def dew_pressure(mixture, thermo, options=None) -> tuple[float, np.ndarray, bool]: ...
def dew_temperature(mixture, thermo, options=None) -> tuple[float, np.ndarray, bool]: ...
```

Algorithm: Wilson correlation initial K + successive substitution with EOS K-values.

### Stability Test (`flash/stability.py`)

```python
def stability_test(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[list[int], float, float, dict]:
    """Returns (stability_flags, SL, SV, result_dict)."""
```

`result_dict` keys: `overall` (`'stable'|'unstable'|'inconclusive'`), `message`.

---

## Step 8: Tools (Convenience API)

`tools.py` — mirrors MATLAB `Tools/` wrappers:

```python
def add_components(names: list[str]) -> list[Component]: ...
def add_mixture(components: list[Component], T: float, p: float) -> Mixture: ...
def add_thermo() -> ThermoModel: ...
def zero_bip(components: list[Component]) -> BIP: ...
```

---

## Step 9: Public API (`src/pvttool/__init__.py`)

```python
from pvttool.classes.component import Component
from pvttool.classes.bip import BIP
from pvttool.classes.mixture import Mixture
from pvttool.classes.thermo_model import ThermoModel
from pvttool.classes.flash_options import FlashOptions
from pvttool.eos.pr import preos
from pvttool.eos.srk import srkeos
from pvttool.eos.pr78 import pr78eos
from pvttool.flash.vle import vle_flash, vle_flash_negative
from pvttool.flash.lle import lle_flash
from pvttool.flash.saturation import bubble_pressure, bubble_temperature, dew_pressure, dew_temperature
from pvttool.flash.stability import stability_test, stability_lle_test
from pvttool.activity.nrtl import nrtl
from pvttool.activity.wilson import wilson
from pvttool.activity.uniquac import uniquac
from pvttool.tools import add_components, add_mixture, add_thermo, zero_bip
```

**Equivalent quick-start (mirrors MATLAB `methanol_water_vle.m`):**

```python
import numpy as np
from pvttool import add_components, add_mixture, add_thermo, FlashOptions, vle_flash

comps = add_components(["Methanol", "Water"])
mix = add_mixture(comps, T=350.0, p=101325.0)
mix.mole_fraction = np.array([0.5, 0.5])
thermo = add_thermo()
opts = FlashOptions()
y, x, V = vle_flash(mix, thermo, opts)
print(f"Vapor fraction: {V:.4f}")
print(f"Vapor composition: {y}")
print(f"Liquid composition: {x}")
```

---

## Step 10: Naming Convention

| MATLAB                     | Python                       |
|----------------------------|------------------------------|
| `vleflash`                 | `vle_flash`                  |
| `vleflashnegative`         | `vle_flash_negative`         |
| `lleflash`                 | `lle_flash`                  |
| `bubblePressure`           | `bubble_pressure`            |
| `bubbleTemperature`        | `bubble_temperature`         |
| `dewPressure`              | `dew_pressure`               |
| `dewTemperature`           | `dew_temperature`            |
| `stabilityTest`            | `stability_test`             |
| `stabilityLLETest`         | `stability_lle_test`         |
| `addComponents`            | `add_components`             |
| `addMixture`               | `add_mixture`                |
| `addThermo`                | `add_thermo`                 |
| `zeroBIP`                  | `zero_bip`                   |
| `Component.fromDatabase`   | `Component.from_database`    |
| `ThermoModel.EOS`          | `ThermoModel.eos`            |
| `ThermoModel.mixingrule`   | `ThermoModel.mixing_rule`    |
| `FlashOptions.maxIteration`| `FlashOptions.max_iteration` |

---

## Step 11: Docstrings

Every public class, method, and function gets a NumPy-style docstring. Example:

```python
def vle_flash(
    mixture: Mixture,
    thermo: ThermoModel,
    options: FlashOptions | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Isothermal vapor-liquid equilibrium flash calculation.

    Solves the Rachford-Rice equation for the vapor fraction and phase
    compositions at fixed temperature and pressure using successive
    substitution with GDEM acceleration.

    Parameters
    ----------
    mixture : Mixture
        Multicomponent mixture with temperature, pressure, and overall
        mole fractions set.
    thermo : ThermoModel
        Thermodynamic model configuration (EOS, mixing rule, etc.).
    options : FlashOptions, optional
        Convergence tolerances and iteration limits. Uses defaults if None.

    Returns
    -------
    vapor_y : np.ndarray, shape (n,)
        Vapor-phase mole fractions.
    liquid_x : np.ndarray, shape (n,)
        Liquid-phase mole fractions.
    vapor_fraction : float
        Molar vapor fraction in [0, 1].
    """
```

---

## Step 12: Tests

Use `pytest`. Reference values come from the existing MATLAB test suite
and examples.

| File                      | MATLAB equivalent         | What is asserted                             |
|---------------------------|---------------------------|----------------------------------------------|
| `test_classes.py`         | `test_classes.m`          | Field types, BIP zero init, from_database    |
| `test_eos.py`             | `test_eos.m`              | Z-factors, residual HR for pure CH₄ at known T/P |
| `test_flash_vle.py`       | `test_flash_vle.m`        | y, x, V for methanol/water; mass balance     |
| `test_flash_lle.py`       | `test_stability_lle.m`    | Two-liquid equilibrium compositions          |
| `test_saturation.py`      | `test_saturation.m`       | Bubble/dew P and T, convergence flags        |
| `test_stability.py`       | `test_stability.m`        | Flags for known stable and unstable systems  |
| `test_residual_props.py`  | `test_residual_props.m`   | HR, SR, GR, VR, Cp_R, Cv_R signs and values |

Tolerances: `atol=1e-5` for mole fractions; `atol=1.0` for pressures/enthalpies.

---

## Step 13: Implementation Order

1. `scripts/convert_puredata.py` + generate `src/pvttool/data/puredata.json`
2. `pyproject.toml` + directory scaffold + `__init__.py` stubs
3. `classes/` — Component, BIP, Mixture, ThermoModel, FlashOptions
4. `eos/_roots.py` — `select_z_roots()`
5. `eos/pr.py`, `eos/srk.py`, `eos/pr78.py`
6. `mixing_rules.py`
7. `activity/nrtl.py`, `activity/wilson.py`, `activity/uniquac.py`
8. `flash/_kvalue.py`, `flash/_rachford_rice.py`, `flash/_fugacity.py`
9. `flash/vle.py`, `flash/lle.py`
10. `flash/saturation.py`, `flash/stability.py`
11. `tools.py`, `auxiliary.py`
12. Populate `src/pvttool/__init__.py` exports
13. `tests/` — one file per module, all passing
14. Docstring pass over all public API

---

## Verification

```bash
# Install in editable mode
uv pip install -e ".[dev]"

# Run tests
uv run pytest tests/ -v --tb=short

# Smoke test
uv run python -c "
import numpy as np
from pvttool import add_components, add_mixture, add_thermo, vle_flash

comps = add_components(['Methanol', 'Water'])
mix = add_mixture(comps, 350.0, 101325.0)
mix.mole_fraction = np.array([0.5, 0.5])
y, x, V = vle_flash(mix, add_thermo())
print(f'V={V:.4f}  y={y}  x={x}')
"
```

Expected: results match MATLAB `methanol_water_vle.m` within `1e-5` mole fraction.
