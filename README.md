# PVTtool

A Python package for PVT (Pressure-Volume-Temperature) calculations using cubic equations of state. Supports vapor-liquid equilibrium (VLE), liquid-liquid equilibrium (LLE), phase stability testing, bubble/dew point calculations, and residual enthalpy calculations for multicomponent mixtures.

## Installation

```bash
uv pip install -e ".[dev]"
```

Requires Python 3.10+ and numpy.

## Quick Start

```python
from pvt import Component, Mixture, ThermoModel, FlashOptions, PREOS
from pvt import vleflash, stabilityTest, vleflashnegative

# 1. Load components from the built-in database
comps, flag = Component.from_database_array(['CH4', 'C2H6', 'C10H22'])

# 2. Create a mixture at T = 300 K, P = 5 MPa
mix = Mixture(comps, 300, 5e6)
mix.mole_fraction = [0.6, 0.2, 0.2]

# 3. Configure thermodynamic model (defaults: PR EOS, vdW mixing)
thermo = ThermoModel()

# 4. Run a stability test
flag, SL, SV, result = stabilityTest(mix, thermo)
print(result["message"])

# 5. Run a VLE flash
opts = FlashOptions()
y, x, V = vleflash(mix, thermo, opts)
print(f"Vapor fraction: {V:.4f}")
```

## Running Tests

```bash
pytest tests/ -v
```

## API Reference

### Classes

| Class | Purpose | Constructor |
|---|---|---|
| `Component` | Pure-component thermodynamic properties | `Component.from_database('CH4')` |
| `BIP` | Binary interaction parameter matrices | `BIP(n)` |
| `Mixture` | Multicomponent mixture at T, P | `Mixture(components, T_K, p_Pa)` |
| `ThermoModel` | EOS + mixing rule + activity model selection | `ThermoModel()` |
| `FlashOptions` | Convergence settings for flash/stability | `FlashOptions()` |

### ThermoModel fields

| Field | Default | Options |
|---|---|---|
| `EOS` | `PREOS` | `PREOS`, `SRKEOS`, `PR78EOS` |
| `activity_model` | `NRTL` | `NRTL`, `Wilson`, `UNIQUAC`, `Margules2` |
| `mixingrule` | `1` | 1=vdW, 2=HV, 3=MHV1, 4=MHV2 |
| `phase` | `1` | 1=liquid Z-root, 2=vapor Z-root |
| `fugacity_switch` | `1` | 1=compute fugacity, 0=skip |

### Flash functions

| Function | Description |
|---|---|
| `vleflash(mix, thermo, opts)` | VLE flash, vapor fraction in [0, 1]; GDEM acceleration |
| `vleflashnegative(mix, thermo, opts)` | VLE flash, extended vapor fraction |
| `lleflash(mix, thermo, opts)` | LLE flash (two liquid phases) |
| `bubblePressure(mix, thermo, opts)` | Bubble-point pressure at fixed T |
| `bubbleTemperature(mix, thermo, opts)` | Bubble-point temperature at fixed P |
| `dewPressure(mix, thermo, opts)` | Dew-point pressure at fixed T |
| `dewTemperature(mix, thermo, opts)` | Dew-point temperature at fixed P |
| `stabilityTest(mix, thermo)` | Michelsen VLE stability test |
| `stabilityLLETest(mix, thermo)` | Michelsen LLE stability test |

### EOS functions

| Function | Description |
|---|---|
| `PREOS(mix, thermo)` | Peng-Robinson (1976) |
| `SRKEOS(mix, thermo)` | Soave-Redlich-Kwong (1972) |
| `PR78EOS(mix, thermo)` | Peng-Robinson with 1978 alpha correction |

Set `compute_props=True` to get residual properties as 5th return value (ResidualProps namedtuple with HR, SR, GR, VR, Cp_R, Cv_R).

### Convenience wrappers

| Function | Description |
|---|---|
| `add_components(names)` | Load Component array from database |
| `add_mixture(comp, T, p)` | Create Mixture |
| `add_thermo()` | Create ThermoModel with defaults |
| `zero_bip(comp)` | Create zero BIP |

## License

BSD 2-Clause — see source file headers.
Copyright 2012-2013 Ali Akbar Eftekhari.
