# PVTtool

A MATLAB/Octave toolbox for PVT (Pressure-Volume-Temperature) calculations using cubic equations of state. Supports vapor-liquid equilibrium (VLE), liquid-liquid equilibrium (LLE), phase stability testing, and residual enthalpy calculations for multicomponent mixtures.

---

## Getting Started

```matlab
% 1. Initialise (adds all subdirectories to path)
PVTinitialize()

% 2. Load components from the built-in database
[comp, flag] = addComponents({'CH4', 'C2H6', 'C10H22'});

% 3. Create a mixture at T = 300 K, P = 5 MPa
mix = addMixture(comp, 300, 5e6);
mix.mole_fraction = [0.6 0.2 0.2];

% 4. Configure thermodynamic model (defaults: PR EOS, vdW mixing)
thermo = addThermo();

% 5. Run a stability test
[flag, SL, SV, result] = stabilityTest(mix, thermo);
disp(result.message)

% 6. Run a VLE flash
opts = FlashOptions();
[y, x, V] = vleflash(mix, thermo, opts);
fprintf('Vapor fraction: %.4f\n', V)
```

---

## Class Reference

All classes are defined in `Classes/` and are available after `PVTinitialize()`.

| Class | Purpose | Constructor |
|---|---|---|
| `Component` | Pure-component thermodynamic properties | `Component(name, formula, ...)` or `Component.fromDatabase('CH4')` |
| `BIP` | Binary interaction parameter matrices | `BIP(n)` — initialises all matrices to zero for n components |
| `Mixture` | Multicomponent mixture at T, P | `Mixture(components, T_K, p_Pa)` |
| `ThermoModel` | EOS + mixing rule + activity model selection | `ThermoModel()` — then override fields |
| `FlashOptions` | Convergence settings for flash/stability | `FlashOptions()` — then override fields |

### ThermoModel fields

| Field | Default | Options |
|---|---|---|
| `EOS` | `@PREOS` | `@PREOS`, `@SRKEOS`, `@PR78EOS` |
| `activity_model` | `@NRTL` | `@NRTL`, `@Wilson`, `@UNIQUAC`, `@Margules2` |
| `mixingrule` | `1` | 1=van der Waals, 2=Huron-Vidal, 3=MHV1, 4=MHV2 |
| `phase` | `1` | 1=liquid Z-root, 2=vapor Z-root |
| `fugacity_switch` | `1` | 1=compute fugacity, 0=skip |

### FlashOptions fields

| Field | Default | Description |
|---|---|---|
| `accuracy` | `1e-7` | Convergence tolerance for VLE/LLE flash |
| `iteration` | `100` | Maximum successive substitution iterations |
| `trivialSolutionMaxError` | `1e-5` | Stability test trivial solution threshold |
| `convergenceMaxError` | `1e-10` | Stability test convergence threshold |
| `maxIteration` | `50` | Maximum stability test iterations |

---

## Setting Binary Interaction Parameters

The `Mixture` constructor creates a fully-initialised `BIP` object with all matrices set to zero. Set non-zero parameters by assigning individual fields:

```matlab
mix = addMixture(comp, T, p);

% EOS kij parameters
mix.bip.EOScons = [0 -0.05; -0.05 0];

% NRTL parameters  (A + B*T + C*T² + D/T + E*ln(T))
mix.bip.NRTLcons    = [0 A12; A21 0];     % [J/mol]
mix.bip.NRTLtdep    = [0 B12; B21 0];
mix.bip.NRTLalfa    = [0 0.3; 0.3 0];     % non-randomness
```

All BIP fields: `EOScons`, `EOStdep`, `NRTLcons`, `NRTLtdep`, `NRTLtdep2`, `NRTLtdepm1`, `NRTLtdeplog`, `NRTLalfa`, `Wilsoncons`, `Wilsontdep`, `UNIQUACcons`, `UNIQUACtdep`, `UNIQUACR`, `UNIQUACQ`.

---

## API Reference

### Flash functions

| Function | Description |
|---|---|
| `vleflash(mix, thermo, opts)` | VLE flash, vapor fraction ∈ [0, 1] |
| `vleflashnegative(mix, thermo, opts)` | VLE flash, extended vapor fraction (negative saturation) |
| `lleflash(mix, thermo, opts)` | LLE flash (two liquid phases) |
| `bubbleTemperature(mix, thermo, opts)` | Bubble-point temperature at fixed P |
| `dewTemperature(mix, thermo, opts)` | Dew-point temperature at fixed P |

### Stability tests

| Function | Description |
|---|---|
| `stabilityTest(mix, thermo)` | Michelsen VLE stability test; returns `(flag, SL, SV, result)` |
| `stabilityLLETest(mix, thermo)` | Michelsen LLE stability test; returns `(flag, SL, SV, result)` |

`stability_flag` values: `1` = stable (trivial), `2` = unstable (non-trivial), `3` = inconclusive.  
`result.overall`: `'stable'` | `'unstable'` | `'inconclusive'`.

### EOS functions

All EOS functions share the signature `[liquid_z, vapor_z, fugacity, HR] = EOS(mixture, thermo)`.

| Function | Description |
|---|---|
| `PREOS` | Peng-Robinson (1976) |
| `SRKEOS` | Soave-Redlich-Kwong (1972) |
| `PR78EOS` | Peng-Robinson with 1978 alpha correction (better for ω > 0.491) |

### Activity models

All share the signature `[gErt, gama] = Model(T, x, component, bip)`.

| Function | Description |
|---|---|
| `NRTL` | Non-Random Two-Liquid |
| `Wilson` | Wilson model |
| `UNIQUAC` | Universal Quasi-Chemical |
| `Margules2` | Two-parameter Margules |

### Mixing rules

| Number | Name | Notes |
|---|---|---|
| 1 | van der Waals | No activity model needed |
| 2 | Huron-Vidal (HV) | Requires activity model |
| 3 | Modified HV-1 (MHV1) | Requires activity model |
| 4 | Modified HV-2 (MHV2) | Requires activity model |

### Convenience functions (Tools/)

| Function | Description |
|---|---|
| `addComponents(names)` | Load Component array from database |
| `addMixture(comp, T, p)` | Create Mixture (wraps `Mixture()`) |
| `addThermo()` | Create ThermoModel with defaults (wraps `ThermoModel()`) |
| `zeroBIP(comp)` | Create zero BIP (wraps `BIP(n)`) |

---

## Examples

All examples are in the `Examples/` folder and can be run after `PVTinitialize()`.

| File | System | Calculation |
|---|---|---|
| `methanol_water_vle.m` | CH₃OH + H₂O at 322.91 K | VLE flash vs experimental data |
| `testcase2.m` | CH₄ to C₁₅H₃₂ (6 components) | Multi-component VLE |
| `testcase3.m` | DME + H₂O + C₁₀H₂₂ | LLE flash + stability test |
| `pure_component_density.m` | H₂O | Liquid/vapor molar volume vs steam tables |
| `bubbledewtest.m` | Binary mixture | Bubble/dew point via fzero |
| `blackoilmodel.m` | Black-oil system | Black-oil thermodynamic model |

---

## Tests

Run the test suite from the project root:

```matlab
PVTinitialize()
cd Tests
run_all_tests
```

Tests cover: class construction, EOS Z-factors, VLE flash vs experimental data, and stability test stable/unstable cases.

---

## Octave Compatibility

Tested with Octave 7+. All classes use value semantics (no `handle`). No MATLAB-only toolboxes are required.

---

## License

BSD 2-Clause — see source file headers.  
Copyright © 2012–2013 Ali Akbar Eftekhari.
