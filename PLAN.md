# PVTtool Python Conversion Plan

## Overview

Convert PVTtool from MATLAB/Octave to a modern Python package using `uv` for environment management, `numpy` for numerical arrays, and `pytest` for testing. The Python API mirrors the MATLAB API while following Python conventions (snake_case, dataclasses, type hints).

## Package structure

```
pvt/
  __init__.py             # Public API + convenience wrappers
  _constants.py           # Physical constants (R = 8.314)
  _data.py                # Component database loader (puredata.json)
  _equations.py           # Vapor-pressure/heat-capacity correlation functions
  component.py            # Component dataclass + from_database()
  bip.py                  # BIP dataclass (all interaction matrices)
  mixture.py              # Mixture dataclass (T, P, mole_fraction, bip)
  thermo_model.py         # ThermoModel dataclass (EOS, activity, mixing, phase)
  flash_options.py        # FlashOptions dataclass (tolerances, max iterations)
  utils.py                # mynormalize(), bisection()
  mixing_rules.py         # mixing_rule() — vdW, HV, MHV1, MHV2 (+ incomplete Wong-Sandler)
  activity.py             # NRTL, Wilson, UNIQUAC, Margules2 activity models
  eos.py                  # select_z_roots(), PREOS(), SRKEOS(), PR78EOS()
  kvalues.py              # kwal_estimate(), kvalue(), kvalueLLE(), fugacity()
  rachford_rice.py        # RachfordRiceNR() — bounded Newton-Raphson
  flash.py                # vleflash(), vleflashnegative(), lleflash(), stabilityTest(), stabilityLLETest()
  saturation.py           # bubbleTemperature(), bubblePressure(), dewTemperature(), dewPressure()
tests/
  test_classes.py
  test_eos.py
  test_flash_vle.py
  test_stability.py
  test_stability_lle.py
  test_residual_props.py
  test_saturation.py
data/
  puredata.json           # Exported from puredata.mat (one-time conversion)
```

## Dependencies

- `numpy` — array operations, roots(), vectorized math
- `pytest` — test runner (dev dependency)

No other dependencies. `scipy` is NOT required (the EOS uses `numpy.roots()` directly).

## Key design decisions

| MATLAB | Python |
|---|---|
| `[1 x N]` row vector | `np.ndarray` shape `(n,)` — all composition/K vectors are 1D |
| `N×N` BIP matrices | `np.ndarray` shape `(n, n)`, symmetric |
| `@PREOS` function handle | Plain function object `PREOS` |
| `nargout >= 5` optional output | Return `(zL, zV, fug, HR)` always; `PREOS_with_props()` returns `PVTResult` NamedTuple |
| `load puredata.mat` | `_data.py` loads from committed `puredata.json` |
| `roots()` for cubic | `numpy.roots()` (same algorithm) |
| `strcmpi` | `str.casefold()` |
| `error('notImplemented')` | `raise NotImplementedError(...)` |
| Value classes (copy on assign) | `@dataclass(frozen=True)` |
| `feval(name)` test runner | `pytest` |
| `[component.Tc]` bracket expansion | `np.array([c.Tc for c in components])` |

## Implementation order

1. Project scaffold (`pyproject.toml`, directory structure)
2. Export `puredata.mat` → `puredata.json` (one-time MATLAB script)
3. `_constants.py`, `_data.py`, `_equations.py`
4. `component.py`, `bip.py`, `mixture.py`, `thermo_model.py`, `flash_options.py`
5. `utils.py`, `mixing_rules.py`, `activity.py`
6. `eos.py` (select_z_roots + PREOS, SRKEOS, PR78EOS)
7. `kvalues.py`, `rachford_rice.py`
8. `flash.py`, `saturation.py`
9. `__init__.py` with public API + convenience wrappers
10. Test suite (pytest)
11. Delete MATLAB files, run full test suite

## Known limitations preserved from MATLAB

- Mixing rule 5 (Wong-Sandler) raises `NotImplementedError`
- Margules2 model is a stub
- UNIQUAC requires `r_vol`/`q_area` properties on Component (added)
- Saturation functions use simple successive substitution (no global convergence)
- PREOS `dadT` uses kij=0 approximation

## Verification

Run the test suite and compare numerical results against MATLAB reference values for all existing test cases.
