# AGENTS.md

## Setup

```bash
uv venv
uv pip install -e ".[dev]"
```

Uses `uv` for environment management. Python 3.10+. Only runtime dependency is `numpy`.

## Run tests

```bash
uv run pytest tests/ -v
```

41 tests covering classes, EOS, flash, stability, residual properties, saturation.

## Architecture

| Module | Contents |
|---|---|
| `pvt/component.py` | `Component` dataclass + `from_database()` / `from_database_array()` |
| `pvt/bip.py` | `BIP` class — all 12 interaction matrices zero-initialised |
| `pvt/mixture.py` | `Mixture` dataclass (components, temperature, pressure) |
| `pvt/thermo_model.py` | `ThermoModel` with EOS/activity/mixing defaults |
| `pvt/flash_options.py` | `FlashOptions` with convergence defaults |
| `pvt/eos.py` | `PREOS`, `SRKEOS`, `PR78EOS`, `select_z_roots`, `ResidualProps` |
| `pvt/mixing_rules.py` | `mixing_rule()` — rules 1-4 (5 is NotImplementedError) |
| `pvt/activity.py` | NRTL, Wilson, UNIQUAC, Margules2 stubs |
| `pvt/flash.py` | `vleflash`, `lleflash`, `stabilityTest`, `stabilityLLETest` |
| `pvt/saturation.py` | bubble/dew temperature and pressure |
| `pvt/kvalues.py` | K-value estimation and update functions |
| `pvt/rachford_rice.py` | `RachfordRiceNR` with bounds enforcement |
| `pvt/utils.py` | `mynormalize`, `bisection` |
| `pvt/_constants.py` | `R = 8.314` |
| `pvt/_data.py` | Component database loader (from `data/puredata.json`) |
| `pvt/_equations.py` | Correlation equations (antoine, DIPPR, polynomial, Aly-Lee) |

## Critical conventions

- **All composition vectors are 1D numpy arrays `(n,)`** — never 2D `(1, n)` or `(n, 1)`.
- **SI units throughout**: temperature [K], pressure [Pa], energy [J/mol], heat capacity [J/(mol·K)], volume [m³/mol].
- **BIP matrices are `n×n`** symmetric numpy arrays. `BIP(n)` zero-initialises every field.
- **EOS signatures**: `(mixture, thermo, compute_props=False)` → `(zL, zV, fugacity, HR, props_or_None)`.
  Set `compute_props=True` for the 5th `ResidualProps` namedtuple.
- **ThermoModel** uses function objects, not handles: `thermo.EOS = PREOS`, `thermo.activity_model = NRTL`.
- **Component lookup**: `Component.from_database_array(['CH4', 'H2O'])` returns `(components, not_found_indices)`.
- GDEM acceleration is embedded in `vleflash` (every 5 iterations). Do not add a second acceleration layer.
- Mixing rule 5 (Wong-Sandler) raises `NotImplementedError`.

## Known limitations

- Activity model stubs (Wilson, UNIQUAC, Margules2) return unity activity coefficients.
- NRTL is fully implemented for N=2; N>2 is vectorized.
- PREOS `dadT` uses kij=0 approximation (consistent with MATLAB).
- Wong-Sandler (rule 5) is deliberately incomplete.
