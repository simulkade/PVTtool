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

---

# Phase 3: Flash Algorithm Improvements and Residual Properties

## Context

After completing the struct→classdef refactoring and adding the test suite, the next improvements address two areas:

1. **Flash algorithm robustness** — The current algorithms work but have several known weaknesses: the cubic EOS root selection is fragile near the critical point, bubble/dew point functions are empty stubs, and the Rachford-Rice solver uses an absolute step tolerance that can fail for near-trivial splits.

2. **Residual property calculations** — The package computes residual enthalpy (HR) in PREOS and PR78EOS, but residual entropy (SR), residual volume (VR), residual Gibbs energy (GR), and heat capacity departures (Cp_R, Cv_R) are not implemented in any EOS. SRKEOS returns HR = 0. These properties are needed for entropy balances, compressor/expander work calculations, and speed-of-sound estimates.

---

## Part A: Flash Algorithm Improvements

### A1. EOS Cubic Root Selection

**Current state:** All three EOS files (`PREOS.m`, `SRKEOS.m`, `PR78EOS.m`) use MATLAB's `roots()` to find the three Z-roots of the cubic, then select:
```matlab
if (sum(imag(z_root)~=0)==0)   % all three roots are real
    liquid_z = min(z_root);
    vapor_z  = max(z_root);
else                            % one real root
    liquid_z = z_root(imag(z_root)==0);
    vapor_z  = liquid_z;
```

**Problems:**
- `imag(z_root)==0` can fail due to floating-point noise (roots() returns tiny non-zero imaginary parts even for the real root of a one-real-root cubic).
- No check that `liquid_z > b*p/(R*T)` — small negative or near-zero Z roots can occur near the critical point.
- The `real()` call in `mixing_rule.m:MHV2` silently discards imaginary parts, masking numerical issues.

**Fix — robust root selection helper (`EOS/select_z_roots.m`):**

```matlab
function [zL, zV] = select_z_roots(z_roots, B_coef)
% Returns liquid (min physical) and vapor (max physical) Z-roots.
% Roots with imag/real ratio > tol or real part <= B_coef are rejected.
tol = 1e-6;
real_z = real(z_roots);
is_real = abs(imag(z_roots)) ./ max(abs(real_z), 1) < tol;
physical = is_real & (real_z > B_coef);
phys_z   = sort(real_z(physical));
if numel(phys_z) >= 2
    zL = phys_z(1);
    zV = phys_z(end);
elseif numel(phys_z) == 1
    zL = phys_z(1);
    zV = phys_z(1);
else
    % Fallback: take real part of root with smallest imaginary part
    [~, idx] = min(abs(imag(z_roots)));
    zL = real(z_roots(idx));
    zV = zL;
end
```

**Integration:** Replace the 6-line root-selection block in PREOS.m, SRKEOS.m, and PR78EOS.m with `[liquid_z, vapor_z] = select_z_roots(z_roots, B_coef)`.

**Files:** Create `EOS/select_z_roots.m`; modify `EOS/PREOS.m`, `EOS/SRKEOS.m`, `EOS/PR78EOS.m` (lines ~93-100 in each).

---

### A2. Rachford-Rice Solver Improvements

**Current state (`Flash/RachfordRiceNR.m`):**
- Convergence criterion: `abs(dvf) < 1e-10` (absolute step size).
- No bounds enforcement: `vapor_frac` can drift outside `(Vmin, Vmax)` where `Vmin = -1/(K_max-1)` and `Vmax = 1/(1-K_min)`.
- No initial bracket selection.

**Problems:**
- For near-trivial splits (K ≈ 1), the denominator `1+V*(K-1)` ≈ 1 and the absolute step `dvf` is already tiny even before a meaningful solution is found.
- No bound prevents negative mole fractions in the composition calculation.

**Fix:**

```matlab
function [vf, RRflag] = RachfordRiceNR(composition, ki, vapor_frac_est)
% Solve Rachford-Rice: sum(z*(K-1)/(1+V*(K-1))) = 0  by Newton-Raphson.
% Bounds V to the physical interval (Vmin, Vmax) defined by K-values.
Kmax = max(ki); Kmin = min(ki);
Vmin = 1/(1 - Kmax);   % lower bound (avoids division by zero)
Vmax = 1/(1 - Kmin);   % upper bound
vapor_frac = min(max(vapor_frac_est, Vmin + 1e-10), Vmax - 1e-10);
RRflag = 0; count = 0;
while count < 100
    count = count + 1;
    nz = composition ~= 0;
    denom = 1 + vapor_frac * (ki(nz) - 1);
    f    = sum(composition(nz) .* (ki(nz)-1) ./ denom);
    dfdv = -sum(composition(nz) .* (ki(nz)-1).^2 ./ denom.^2);
    dvf  = f / dfdv;
    vapor_frac = vapor_frac - dvf;
    vapor_frac = min(max(vapor_frac, Vmin + 1e-10), Vmax - 1e-10);
    if abs(dvf) < 1e-10 * (1 + abs(vapor_frac))   % relative tolerance
        RRflag = 1; break
    end
end
vf = vapor_frac;
```

**File:** `Flash/RachfordRiceNR.m` — full replacement of the while-loop body (~lines 46-70).

---

### A3. Bubble and Dew Point Functions

**Current state:** All four functions (`bubblePressure.m`, `bubbleTemperature.m`, `dewPressure.m`, `dewTemperature.m`) are empty stubs.

**Algorithm (successive-substitution for saturation):**

The bubble/dew point conditions are:
- **Bubble pressure** (T fixed, V=0): `Σ K_i * z_i = 1` where `K_i = φ_i^L / φ_i^V`
- **Dew pressure** (T fixed, V=1): `Σ z_i / K_i = 1`
- **Bubble temperature** (P fixed, V=0): same constraint, iterate on T
- **Dew temperature** (P fixed, V=1): same constraint, iterate on T

**Implementation for `bubbleTemperature(mix, thermo, opts)`:**

```
Signature: [T_bub, y_bub, flag] = bubbleTemperature(mix, thermo, opts)

1. Initial K from Wilson correlation (kval_estimate)
2. Initial T from mix.temperature (used as starting guess)
3. Successive substitution loop (max opts.iteration):
   a. Normalise incipient vapor: y = K.*z / sum(K.*z)
   b. Build vapor mixture at current T, p; call thermo.EOS for φ^V
   c. Build liquid mixture at current T, p; call thermo.EOS for φ^L
   d. Update K_new = φ^L ./ φ^V
   e. Check S = sum(K.*z) — at bubble: S → 1
   f. Adjust T via: T = T * (log(S) / sum(z.*log(K)) or bisection on (S-1))
   g. Convergence: |S-1| < opts.accuracy AND max|K_new-K|/K < opts.accuracy
4. Return T, y (normalized incipient vapor), flag (1=converged, 0=not)
```

The pressure analogue (`bubblePressure`) iterates on P instead of T; the dew variants flip liquid↔vapor and use `Σ z/K = 1`.

A clean, self-contained implementation of all four functions shares a helper:

```matlab
% Flash/saturation_ss.m (internal helper, not public API)
function [sat_val, incipient_x, flag] = saturation_ss(mix, thermo, opts, mode)
% mode: 'bub_T','bub_P','dew_T','dew_P'
```

**Files to create:**
- `Flash/bubblePressure.m` — full implementation
- `Flash/bubbleTemperature.m` — full implementation
- `Flash/dewPressure.m` — full implementation
- `Flash/dewTemperature.m` — full implementation

**Signatures (consistent with existing package style):**
```matlab
[T_bub, y, flag] = bubbleTemperature(mix, thermo, opts)
[P_bub, y, flag] = bubblePressure(mix, thermo, opts)
[T_dew, x, flag] = dewTemperature(mix, thermo, opts)
[P_dew, x, flag] = dewPressure(mix, thermo, opts)
```

---

### A4. Flash Convergence Acceleration (GDEM)

**Current state:** `vleflash.m` and `vleflashnegative.m` use plain successive substitution (SS). SS can converge very slowly (hundreds of iterations) near the critical point or for highly non-ideal systems.

**Proposed improvement:** Add Michelsen's GDEM (Dominant Eigenvalue Method) acceleration as an optional step every `n_acc = 5` iterations.

GDEM uses the last three successive-substitution iterates to estimate the dominant eigenvalue `λ` of the Jacobian and jump to the fixed point:

```matlab
% After 3 SS steps, collect log(K) vectors:
% g0 = log(K_n), g1 = log(K_{n-1}), g2 = log(K_{n-2})
dg1 = g1 - g0; dg2 = g2 - g1;
lambda = (dg2' * dg1) / (dg1' * dg1);  % dominant eigenvalue estimate
if 0 < lambda && lambda < 1
    g_acc = g0 - lambda/(1-lambda) * dg1;
    K_acc = exp(g_acc);
end
```

This is a targeted, well-understood acceleration that is compatible with the existing SS loop structure. It should be applied only when `lambda` is in `(0,1)` (contracting iteration), and fall back to plain SS otherwise.

**Files to modify:** `Flash/vleflash.m` (~lines where K is updated in the SS loop), `Flash/vleflashnegative.m` (same pattern).

---

## Part B: Residual Property Calculations

### B1. Derivation Summary (PR EOS)

For the Peng-Robinson EOS with van der Waals mixing rule, the residual properties are:

**Parameters:**
```
a(T) = Σ_i Σ_j x_i x_j (a_i a_j)^0.5 (1-k_ij)
b    = Σ_i x_i b_i
A    = aP/(RT)²     B = bP/(RT)
da/dT = Σ_i Σ_j x_i x_j (a_i a_j)^0.5 (1-k_ij) * (da_i/dT a_j + a_i da_j/dT)/(2(a_i a_j)^0.5 ... )
      = -a * Σ_i x_i m_i/(√(Tc_i) α_i^0.5) / (a^0.5 √T)   [simplified form]
```

The `da/dT` term is already partially computed in the HR calculation in PREOS.m (the `Abar`, `part1` sum). It should be extracted into a shared sub-expression.

**Residual properties (PR, vdW mixing):**

| Property | Formula |
|---|---|
| HR | `RT(Z-1) + (T da/dT - a) / (b√8) * ln((Z+(1+√2)B)/(Z+(1-√2)B))` |
| SR | `R ln(Z-B) + (da/dT) / (b√8) * ln((Z+(1+√2)B)/(Z+(1-√2)B))` |
| GR | `HR - T·SR` (or directly: `RT(Z-1) - RT ln(Z-B) - a/(b√8)*ln(...)`) |
| VR | `RT(Z-1)/P` |
| AR | `GR - P·VR + RT` (Helmholtz, optional) |

**Heat capacities:**

For Cp_R and Cv_R the second temperature derivative of `a` is needed:
```
d²a/dT² = Σ_i x_i * (a_i aci_i) * m_i/(2Tc_i) * (m_i+1) / (T * α_i)
```

Then:
```
Cv_R = -T * d²a/dT² / (b√8) * ln((Z+(1+√2)B)/(Z+(1-√2)B))
Cp_R = Cv_R - R + (RT + (RT² d²a/dT²)(Z-B)/P) ... [full expression from thermodynamic identity]
```

A simpler and numerically stable route for Cp_R uses the thermodynamic identity:
```
Cp_R = Cv_R - R + T*(dP/dV|T)^{-1} * (dP/dT|V)^2 / (P/RT)
```
where `dP/dV` and `dP/dT` can be computed analytically from the cubic EOS.

### B2. Implementation Plan

**Approach:** Add residual properties as additional return values to each EOS function via a `properties` struct, keeping backward compatibility.

**New EOS signature:**
```matlab
[liquid_z, vapor_z, fugacity, HR, props] = PREOS(mixture, thermo)
```

Where `props` is a struct with fields:
- `props.HR`  — residual enthalpy [J/mol] (same as current 4th output)
- `props.SR`  — residual entropy [J/(mol·K)]
- `props.GR`  — residual Gibbs energy [J/mol]
- `props.VR`  — residual molar volume [m³/mol]
- `props.Cp_R` — residual isobaric heat capacity [J/(mol·K)]
- `props.Cv_R` — residual isochoric heat capacity [J/(mol·K)]

Backward compatibility: the 4th output `HR` stays as-is; `props` is only computed/returned when `nargout >= 5`.

**Extraction of `da/dT` as a shared expression:**

Both HR (already computed) and SR require `da/dT`. The current PREOS.m computes this inline as `part1`. This should be extracted before the `if fug_need` block as:

```matlab
% da/dT for vdW mixing (always needed for HR/SR)
daidT = -aci .* mi .* sqrt(Tr) ./ (critical_temp .* alfai);  % da_i/dT
dadT  = 0;
for i = 1:N
    for j = 1:N
        dadT = dadT + x(i)*x(j)*(1-BIP.EOScons(i,j)-BIP.EOStdep(i,j)*T) ...
               * (daidT(i)*sqrt(ai(j)) + sqrt(ai(i))*daidT(j)) / (2*sqrt(ai(i)*ai(j)));
    end
end
% d²a/dT² for Cp_R / Cv_R
d2aidT2 = aci .* mi .* (mi+1) .* Tr ./ (2 * critical_temp.^2 .* alfai);
d2adT2  = 0;  % simplified: only diagonal cross-terms for vdW
for i = 1:N
    d2adT2 = d2adT2 + x(i)^2 * d2aidT2(i);
    % off-diagonal cross-terms omitted for brevity; add if needed
end
```

**SR formula (PR):**
```matlab
ln_term = log((zz + (1+sqrt(2))*B_coef) / (zz + (1-sqrt(2))*B_coef));
SR = R*log(zz - B_coef) + dadT / (b * 2*sqrt(2)) * ln_term;
```

**GR formula:**
```matlab
GR = HR - T * SR;
```

**VR formula:**
```matlab
VR = R*T*(zz - 1) / p;
```

**Cv_R formula:**
```matlab
Cv_R = -T * d2adT2 / (b * 2*sqrt(2)) * ln_term;
```

**Cp_R formula** (via thermodynamic identity, avoids need for full cubic inversion):
```matlab
% dP/dT|V and dP/dV|T at molar volume V = zz*R*T/p
V_mol = zz * R * T / p;
dPdT = R/(V_mol - b) - dadT / (V_mol^2 + 2*b*V_mol - b^2);
dPdV = -R*T/(V_mol-b)^2 + 2*a*(V_mol+b) / (V_mol^2+2*b*V_mol-b^2)^2;
Cp_R = Cv_R - R - T * dPdT^2 / dPdV;
```

### B3. SRK Residual Properties

The SRK EOS uses a different cubic form (Redlich-Kwong-Soave). The analog formulas with `u=1, w=0` (van der Waals-Soave notation) are:

```
HR_SRK = RT(Z-1) + (T da/dT - a)/b * ln(Z/(Z+B))
SR_SRK = R ln(Z-B) + dadT/b * ln(Z/(Z+B))
```

The `da/dT` derivation is the same pattern as PR. SRKEOS.m currently has the HR calculation commented out entirely; it should be implemented using the same structure.

**File:** `EOS/SRKEOS.m` — implement `da/dT`, `HR`, `SR`, `GR`, `VR`, `Cp_R`, `Cv_R`.

---

## Implementation Sequence (Phase 3)

### Priority 1: Residual Properties (high value, self-contained)
1. Add `select_z_roots.m` helper — unblocks numerical stability improvements
2. Modify PREOS.m: extract `da/dT`, implement SR/GR/VR/Cp_R/Cv_R, add optional `props` 5th output
3. Modify PR78EOS.m: same changes (alpha function differs only in `mi`, residual property formulas identical)
4. Modify SRKEOS.m: implement HR (currently 0), SR, GR, VR, Cp_R, Cv_R using SRK analog formulas
5. Add `test_residual_props.m` to the test suite

### Priority 2: Bubble/Dew Points (fills empty stubs)
6. Implement `bubbleTemperature.m` and `bubblePressure.m` with SS algorithm
7. Implement `dewTemperature.m` and `dewPressure.m`
8. Add `test_saturation.m` — test against known bubble/dew temperatures for methanol-water

### Priority 3: Flash Algorithm Robustness
9. Improve `RachfordRiceNR.m`: add bounds, relative tolerance
10. Add GDEM acceleration to `vleflash.m` and `vleflashnegative.m`
11. Update tests to confirm no numerical regressions

---

## Files to Create

| File | Description |
|---|---|
| `EOS/select_z_roots.m` | Robust Z-root selection helper (replaces inline logic) |
| `Tests/test_residual_props.m` | Tests for SR, GR, VR, Cp_R, Cv_R |
| `Tests/test_saturation.m` | Tests for bubble/dew point functions |

## Files to Modify

| File | Change |
|---|---|
| `EOS/PREOS.m` | Extract da/dT; add SR/GR/VR/Cp_R/Cv_R; use select_z_roots; add props 5th output |
| `EOS/PR78EOS.m` | Same as PREOS.m (different mi formula, same residual formulas) |
| `EOS/SRKEOS.m` | Implement HR (currently 0); add SR/GR/VR/Cp_R/Cv_R; use select_z_roots |
| `Flash/RachfordRiceNR.m` | Add bounds enforcement; change to relative+absolute convergence criterion |
| `Flash/vleflash.m` | Add GDEM acceleration every 5 SS steps |
| `Flash/vleflashnegative.m` | Add GDEM acceleration every 5 SS steps |
| `Flash/bubblePressure.m` | Full implementation (currently empty stub) |
| `Flash/bubbleTemperature.m` | Full implementation (currently empty stub) |
| `Flash/dewPressure.m` | Full implementation (currently empty stub) |
| `Flash/dewTemperature.m` | Full implementation (currently empty stub) |

## Verification

```matlab
PVTinitialize()
cd Tests
run_all_tests   % all existing 5 tests must still pass

% Residual properties smoke test
[comp, ~] = addComponents({'CH4'});
mix = Mixture(comp, 300, 5e6);
th = ThermoModel(); th.fugacity_switch = 0;
[zl, zv, ~, HR, props] = PREOS(mix, th);
assert(HR < 0,         'HR negative for liquid CH4')
assert(props.SR < 0,   'SR negative for liquid CH4 (more ordered than ideal gas)')
assert(props.VR < 0,   'VR negative for liquid (higher density than ideal gas)')
assert(props.GR < 0,   'GR negative (stable phase)')

% Bubble temperature test (methanol-water at known conditions)
[comp2, ~] = addComponents({'CH4O','H2O'});
mix2 = Mixture(comp2, 322.91, 26131);
mix2.bip.EOScons = [0 -0.07; -0.07 0];
th2 = ThermoModel();
opts = FlashOptions();
[T_bub, y_bub, flag] = bubbleTemperature(mix2, th2, opts);
assert(flag == 1,                 'bubbleTemperature converged')
assert(abs(T_bub - 322.91) < 5,  'bubble T within 5 K of known value')
```
