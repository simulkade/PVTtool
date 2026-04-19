# PVTtool Refactoring Plan: Struct → classdef

## Context

PVTtool is a MATLAB/Octave package for PVT flash calculations. The core data is currently passed around as plain structs (`mixture`, `thermo_model`, `BIP`, `options`). A `Component` classdef exists but is never actually used — `addComponents` still builds struct arrays. The goal is to replace all structs with proper classdef classes, fix existing bugs discovered in the process, and improve the code where natural to do so, while keeping full Octave compatibility.

---

## Current Struct Inventory

| Struct variable | Created by | Fields |
|---|---|---|
| `component(n)` | `addComponents` | name, formula, MW, Tc, Pc, Vc, Zc, acentric_factor, Psat_eq/coefs, dh_vap_eq/coefs, cp_liq/ig fields, dhf_ig, dgf_ig, ds_ig, dh_comb |
| `mixture` | `addMixture` | components, mole_fraction, pressure, temperature, bip |
| `mixture.bip` | `zeroBIP` | EOScons, EOStdep (only — others missing, see Bug 3) |
| `thermo_model` | `addThermo` | EOS, activity_model, mixingrule, phase, fugacity_switch |
| `options` | user-created | accuracy, iteration (VLE); trivialSolutionMaxError, convergenceMaxError, maxIteration (stability) |

---

## Bugs to Fix During Refactoring

1. **`Component.m:67`** — `obj.name = name` but `name` is not in the constructor parameter list (was likely removed accidentally). The class was never called in practice so this was silent.
2. **`Component.m:125`** — `cp_ig` method calls `obj.cp_ig(T, obj.cp_ig_coefs)` — recursive self-call. Should be `obj.cp_ig_eq(T, obj.cp_ig_coefs)`.
3. **`zeroBIP.m`** — truncated at line 50; only initialises `EOScons` and `EOStdep`. The NRTL model (`NRTL.m:51`) uses `NRTLcons`, `NRTLtdep`, `NRTLtdep2`, `NRTLtdepm1`, `NRTLtdeplog`, `NRTLalfa`. Wilson, UNIQUAC also have BIP fields. Users currently set these manually; the `BIP` class constructor will initialise all fields to zero.
4. **`lleflash.m`** — file content is identical to `vleflash.m` (calls `kvalueLLE` on last iteration but body is otherwise VLE). Needs its own proper LLE implementation using `kvalueLLE`.
5. **`mixing_rule.m:95`** — duplicate `elseif (mixing_rule_num == 4)` block (should be `== 5` for Wong-Sandler); also references undefined variable `BIPs` and `Cstar`. Mark as TODO/incomplete.
6. **`PREOS.m:126`** — Huron-Vidal branch calls `activityfun(T, x, component, BIP)` passing a bare `component` variable that is not in scope (should be `mixture.components`).

---

## New Class Hierarchy

All classes live under `Classes/` in `@ClassName/ClassName.m` format (already on path from `PVTinitialize`).

### 1. `Component` — refactor existing (`Classes/@Component/Component.m`)

**Changes:**
- Add `name` as first constructor parameter (fix Bug 1).
- Fix recursive `cp_ig` call (fix Bug 2).
- Add static method `fromDatabase(name_or_formula)` that loads `puredata.mat` and returns a single `Component` object. This replaces the data-loading logic in `addComponents`.
- Keep all existing properties unchanged.

```matlab
classdef Component
    properties
        name, formula, MW, Tc, Pc, Vc, Zc, acentric_factor
        Psat_eq, Psat_coefs, PsatTrange, Psatrange
        dh_vap_eq, dh_vap_coefs, dh_vap_Trange, dh_vap_range
        cp_liq_eq, cp_liq_coefs, cp_liq_Trange, cp_liq_range
        cp_ig_eq,  cp_ig_coefs,  cp_ig_Trange,  cp_ig_range
        dhf_ig, dgf_ig, ds_ig, dh_comb
    end
    methods
        function obj = Component(name, formula, MW, Tc, Pc, ...) % full list
        function p_sat = vapor_pressure(obj, T)
        function dh_v  = dh_vap(obj, T)
        function cp_l  = cp_liq(obj, T)
        function cp    = cp_ig(obj, T)   % fix: call obj.cp_ig_eq
    end
    methods (Static)
        function comp = fromDatabase(name_or_formula)
        function comps = fromDatabaseArray(names_cell)  % replaces addComponents logic
    end
end
```

### 2. `BIP` — new class (`Classes/@BIP/BIP.m`)

Replaces `zeroBIP`. Constructor initialises all binary interaction parameter matrices to zero for an `n`-component mixture.

```matlab
classdef BIP
    properties
        EOScons, EOStdep
        NRTLcons, NRTLtdep, NRTLtdep2, NRTLtdepm1, NRTLtdeplog, NRTLalfa
        Wilsoncons, Wilsontdep
        UNIQUACcons, UNIQUACtdep, UNIQUACR, UNIQUACQ
    end
    methods
        function obj = BIP(n)
            obj.EOScons  = zeros(n); obj.EOStdep  = zeros(n);
            obj.NRTLcons = zeros(n); obj.NRTLtdep = zeros(n);
            obj.NRTLtdep2 = zeros(n); obj.NRTLtdepm1 = zeros(n);
            obj.NRTLtdeplog = zeros(n); obj.NRTLalfa = zeros(n);
            obj.Wilsoncons = zeros(n); obj.Wilsontdep = zeros(n);
            obj.UNIQUACcons = zeros(n); obj.UNIQUACtdep = zeros(n);
            obj.UNIQUACR = zeros(1,n); obj.UNIQUACQ = zeros(1,n);
        end
    end
end
```

### 3. `Mixture` — refactor existing empty class (`Classes/@Mixture/Mixture.m`)

Replaces `addMixture`. Constructor creates equimolar composition and a `BIP` object.

```matlab
classdef Mixture
    properties
        components       % 1xN Component array
        mole_fraction    % 1xN double
        pressure         % scalar [Pa]
        temperature      % scalar [K]
        bip              % BIP object
    end
    methods
        function obj = Mixture(components, T_K, p_Pa)
            obj.components    = components;
            n                 = numel(components);
            obj.mole_fraction = ones(1,n)/n;
            obj.pressure      = p_Pa;
            obj.temperature   = T_K;
            obj.bip           = BIP(n);
        end
    end
end
```

### 4. `ThermoModel` — new class (`Classes/@ThermoModel/ThermoModel.m`)

Replaces `addThermo`. Stores the thermodynamic model configuration.

```matlab
classdef ThermoModel
    properties
        EOS              = @PREOS
        activity_model   = @NRTL
        mixingrule       = 1    % 1=vdW, 2=HV, 3=MHV1, 4=MHV2, 5=WS
        phase            = 1    % 1=liquid, 2=vapor
        fugacity_switch  = 1
    end
    methods
        function obj = ThermoModel()   % default constructor; user overrides fields
        end
    end
end
```

### 5. `FlashOptions` — new class (`Classes/@FlashOptions/FlashOptions.m`)

Replaces the ad-hoc `options` struct used in flash and stability functions. Provides defaults for both VLE flash and stability test options.

```matlab
classdef FlashOptions
    properties
        accuracy                 = 1e-7
        iteration                = 100
        trivialSolutionMaxError  = 1e-5
        convergenceMaxError      = 1e-10
        maxIteration             = 50
    end
    methods
        function obj = FlashOptions()
        end
    end
end
```

---

## File-by-File Changes

### `Tools/` — Keep as thin backward-compatible wrappers

| File | Change |
|---|---|
| `addComponents.m` | Delegate to `Component.fromDatabaseArray(names)`. Return `Component` array. |
| `addMixture.m` | Delegate to `Mixture(components, T, p)`. |
| `addThermo.m` | Delegate to `ThermoModel()`. |
| `zeroBIP.m` | Delegate to `BIP(n)`. |

These keep existing example scripts working without modification.

### `EOS/PREOS.m`, `EOS/SRKEOS.m`, `EOS/PR78EOS.m`

- Field-access syntax (`mixture.components`, `mixture.bip`, etc.) is identical for structs and objects — no change needed for data extraction lines.
- Fix Bug 6: Huron-Vidal and activity-model branches pass `component` (undefined) — change to `mixture.components`.
- Apply same fix pattern to `SRKEOS.m` and `PR78EOS.m`.

### `MixingRules/mixing_rule.m`

- No signature change needed; `mixture.bip` and `mixture.components` access syntax unchanged.
- Add comment marking Wong-Sandler (rule 5) as unfinished (fix Bug 5).

### `ActivityModels/NRTL.m`, `Wilson.m`, `UNIQUAC.m`, `Margules2.m`

- Signatures stay `(temperature, x, component, BIP)`.
- Property access `[BIP.NRTLcons]` etc. on an object works identically to a struct — no code change required.
- NRTL already handles `NRTLtdep2`, `NRTLtdepm1`, `NRTLtdeplog`; with `BIP` class these are always initialised, so no nil-checks needed.

### `Flash/` functions

- All accept `(mixture, thermo, ...)` — no signature changes.
- `.mole_fraction`, `.pressure`, `.temperature`, `.components` access is unchanged.
- `options.accuracy`, `options.iteration` etc. unchanged (FlashOptions exposes same property names).
- Fix `lleflash.m` (Bug 4) so it uses a proper two-liquid iteration via `kvalueLLE` throughout.

### `PVTinitialize.m`

- Already adds `Classes/` to path. Verify it recurses into `@ClassName` subdirectories — in MATLAB/Octave, adding the parent directory is sufficient; the `@` folders are found automatically.

---

## Octave Compatibility Notes

- Use **value classes** only (no `< handle`) — fully supported in Octave 7+.
- Do not use `events`, `listeners`, `meta.*`, or `notify` — not needed here.
- `methods (Static)` is supported in Octave.
- Array indexing `component_array(n).field` and bracket expansion `[component_array.Tc]` work with classdef object arrays in both MATLAB and Octave.
- Default property values in `properties` blocks (e.g. `EOS = @PREOS`) work in both; if Octave has issues, move defaults into the constructor body.

---

## Migration Order (Recommended Implementation Sequence)

1. **`BIP` class** — no dependencies, unblocks Mixture
2. **`Component` class** — fix bugs, add static `fromDatabase`
3. **`Mixture` class** — depends on Component + BIP
4. **`ThermoModel` class**
5. **`FlashOptions` class**
6. **Update `Tools/`** wrappers to delegate to new classes
7. **Fix EOS bugs** (Bug 6 in PREOS/SRKEOS/PR78EOS)
8. **Fix `zeroBIP`** wrapper (now complete with all fields)
9. **Fix `lleflash.m`** (Bug 4)
10. **Update `PVTinitialize.m`** if needed

---

## Verification

Run all example scripts in `Examples/` before and after refactoring. Key tests:

```matlab
PVTinitialize()

% 1. Smoke test: component loading
[component, flag] = addComponents({'CH4', 'H2O'})
assert(isa(component(1), 'Component'))

% 2. VLE flash (methanol/water)
run('Examples/methanol_water_vle.m')   % produces a plot

% 3. Multi-component VLE
run('Examples/testcase2.m')

% 4. LLE flash
run('Examples/testcase3.m')

% 5. Bubble/dew point
run('Examples/bubbledewtest.m')

% 6. Stability test
mixture1 = Mixture(component, 322.91, 15932);
thermo1  = ThermoModel();
[stability_flag, SL, SV] = stabilityTest(mixture1, thermo1)
```

Expected: all results numerically identical to pre-refactor outputs.

---

## Out of Scope

- EOS and activity model functions remain as standalone functions (not classes). They are used as function handles (`@PREOS`, `@NRTL`) which is clean and already works well.
- Transport properties (`puregasviscosity`, `pureliquidviscosity`) — not part of the struct refactoring, leave as-is.
- The `puredata.mat` database — no changes.
- Completing stub/unfinished functions (`bubblePressure`, `dewPressure`, Wong-Sandler rule) — out of scope for this refactoring.

---

## Files to Create

| File | Type |
|---|---|
| `Classes/@BIP/BIP.m` | New class |
| `Classes/@ThermoModel/ThermoModel.m` | New class |
| `Classes/@FlashOptions/FlashOptions.m` | New class |
| `PLAN.md` | This document (copy to project root) |

## Files to Modify

| File | Change summary |
|---|---|
| `Classes/@Component/Component.m` | Fix constructor (name param), fix cp_ig bug, add static fromDatabase methods |
| `Classes/@Mixture/Mixture.m` | Full implementation (was empty) |
| `Tools/addComponents.m` | Delegate to Component.fromDatabaseArray |
| `Tools/addMixture.m` | Delegate to Mixture constructor |
| `Tools/addThermo.m` | Delegate to ThermoModel constructor |
| `Tools/zeroBIP.m` | Delegate to BIP constructor |
| `EOS/PREOS.m` | Fix component variable scope in activity model branches |
| `EOS/SRKEOS.m` | Same fix as PREOS |
| `EOS/PR78EOS.m` | Same fix as PREOS |
| `Flash/lleflash.m` | Fix to be a proper LLE implementation |
| `MixingRules/mixing_rule.m` | Clarify/comment Wong-Sandler incomplete block |

---

# Phase 2: Tests, Documentation, and Stability Improvements

## Context

After the struct→classdef refactoring, the package needs:
1. A proper test suite (none exists) to guard against regressions.
2. Improved docstrings on the public API so `help` gives useful output.
3. An expanded README covering the new class-based API.
4. `stabilityTest` and `stabilityLLETest` improved to report results in human-readable form in addition to numeric flags, and to fix a logic bug in the post-loop analysis.

---

## 1. Stability Test Improvements

### Bug: post-loop flag analysis is wrong

Both `stabilityTest.m` and `stabilityLLETest.m` contain this pattern:

```matlab
% WRONG (current):
if (conv_error > convergence_eps)    % flag = 1
elseif (triv_error > trivial_eps)    % flag = 2
elseif (j >= max_itr)               % flag = 3
end
```

The while-loop condition is `conv_error>eps AND triv_error>eps AND j<max_itr`.  
It exits when any one condition fails, so **the first `if` above conflates trivial-solution exit with max-iteration exit**.

Correct logic:

```matlab
% CORRECT:
if (triv_error <= trivial_eps)          % exited because trivial → stable
    flag = 1;
elseif (conv_error <= convergence_eps)  % converged to non-trivial → unstable
    flag = 2;
else                                    % j >= max_itr, no convergence
    flag = 3;
end
```

Flag semantics (both VLE and LLE):
| Flag | Meaning |
|---|---|
| 1 | Trivial solution — mixture is stable in this direction |
| 2 | Non-trivial convergence — mixture is unstable (split will occur) |
| 3 | Maximum iterations reached — result is inconclusive |

### Add `result` as optional 4th return value

New signature (backward compatible):

```matlab
function [stability_flag, SL, SV, result] = stabilityTest(mixture, thermo, varargin)
```

`result` is a struct:
```
result.overall         — 'stable' | 'unstable' | 'inconclusive'
result.message         — human-readable summary, e.g.:
                         "Mixture is stable. No phase split detected."
                         "Mixture is UNSTABLE. Vapor-liquid split expected."
result.test1_message   — per-test string
result.test2_message   — per-test string
result.SL              — liquid saturation (same as SL output)
result.SV              — vapor saturation (same as SV output)
```

Overall determination rule (same for VLE and LLE):
- Any test returns flag 2 → `overall = 'unstable'`
- Both tests return flag 1 → `overall = 'stable'`
- No flag 2 but a flag 3 present → `overall = 'inconclusive'`

When the function is called with `nargout == 0` (interactive use without output capture), print `result.message` to the console automatically.

**Files to modify:** `Flash/stabilityTest.m`, `Flash/stabilityLLETest.m`

---

## 2. Test Suite

### Location: `Tests/`

All tests use `assert()` (available in both MATLAB and Octave) via a thin local helper that catches assertion failures and counts results. The runner prints a coloured (or plain) PASS/FAIL summary.

### Test helper: `Tests/pvt_assert.m`

```matlab
function pvt_assert(condition, test_name)
% Throws a named error if condition is false.
if ~condition
    error('PVTtest:failed', 'FAIL: %s', test_name);
end
```

### Runner: `Tests/run_all_tests.m`

Calls each test file inside a `try/catch`, accumulates pass/fail counts, prints summary. No external toolbox required.

### Test files and coverage

| File | What is tested | Reference values |
|---|---|---|
| `test_classes.m` | BIP, Component, Mixture, ThermoModel, FlashOptions constructors; field values; static Component.fromDatabase | Structural checks only |
| `test_eos.m` | PREOS / SRKEOS / PR78EOS: Z-factor, fugacity | Pure CH4 at 300 K, 100 bar: known PR Z-liquid and Z-vapor; cross-check PREOS vs PR78EOS for ω<0.491 (should match) |
| `test_flash_vle.m` | vleflash, vleflashnegative: VLE compositions | Methanol-water at 322.91 K, first three pressure points from `methanol_water_vle.m` experimental data; tolerance 0.01 mole fraction |
| `test_stability.m` | stabilityTest (VLE): stable and unstable cases; new `result` struct fields | Pure CH4 at 300 K/1 bar → stable (both flags = 1); CH4/C10H22 50:50 at 300 K/5 MPa → unstable (a flag = 2) |
| `test_stability_lle.m` | stabilityLLETest (LLE) | Water/decane at 300 K/100 bar → unstable (immiscible) |

---

## 3. Documentation

### 3a. README.md — full rewrite

Structure:
1. **Overview** — what PVTtool does, EOS family supported
2. **Getting Started** — `PVTinitialize`, quick-start code snippet
3. **Class Reference** — table: class, purpose, constructor
4. **API Reference** — tables for flash functions, stability functions, EOS functions, activity models
5. **Mixing Rules** — numbered list (1–4) with names
6. **Examples** — brief description of each file in `Examples/`
7. **Octave Compatibility** — version note
8. **License**

### 3b. Function docstrings

Update SYNOPSIS/PARAMETERS/RETURNS blocks in:
- `Flash/vleflash.m`
- `Flash/vleflashnegative.m`
- `Flash/lleflash.m`
- `Flash/stabilityTest.m`
- `Flash/stabilityLLETest.m`
- `EOS/PREOS.m`
- `EOS/SRKEOS.m`
- `EOS/PR78EOS.m`

---

## Implementation Sequence

1. Fix stability test logic bug + add `result` output in both `stabilityTest.m` and `stabilityLLETest.m`
2. Write `Tests/pvt_assert.m` and `Tests/run_all_tests.m`
3. Write test files (`test_classes.m`, `test_eos.m`, `test_flash_vle.m`, `test_stability.m`, `test_stability_lle.m`)
4. Update function docstrings in Flash/ and EOS/
5. Rewrite README.md

## New Files

| File | Type |
|---|---|
| `Tests/run_all_tests.m` | Test runner |
| `Tests/pvt_assert.m` | Test helper |
| `Tests/test_classes.m` | Unit tests for classdef objects |
| `Tests/test_eos.m` | Unit tests for EOS functions |
| `Tests/test_flash_vle.m` | Integration tests for VLE flash |
| `Tests/test_stability.m` | Tests for VLE stability test |
| `Tests/test_stability_lle.m` | Tests for LLE stability test |

## Modified Files

| File | Change |
|---|---|
| `Flash/stabilityTest.m` | Fix flag logic; add `result` 4th output; add console print |
| `Flash/stabilityLLETest.m` | Same fixes as stabilityTest |
| `Flash/vleflash.m` | Improve docstring |
| `Flash/vleflashnegative.m` | Improve docstring |
| `Flash/lleflash.m` | Improve docstring |
| `EOS/PREOS.m` | Improve docstring |
| `EOS/SRKEOS.m` | Improve docstring |
| `EOS/PR78EOS.m` | Improve docstring |
| `README.md` | Full rewrite |
