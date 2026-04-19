function test_stability_lle()
% test_stability_lle  Tests for stabilityLLETest (LLE).
%
% Cases:
%   1. Water/decane 50:50 at 300 K / 100 bar → immiscible → UNSTABLE.
%   2. result struct is fully populated.
%   3. stabilityLLETest accepts FlashOptions.

thermo = ThermoModel();

% ---- Case 1: H2O / C10H22 — fully immiscible → UNSTABLE -----------------
[comp, ~] = addComponents({'H2O', 'C10H22'});
mix = Mixture(comp, 300, 100e5);
mix.mole_fraction = [0.5 0.5];

[flag, SL, SV, res] = stabilityLLETest(mix, thermo);

pvt_assert(numel(flag) == 2,           'stabilityLLETest: flag has 2 elements');
pvt_assert(SL >= 0,                    'stabilityLLETest: SL >= 0');
pvt_assert(SV >= 0,                    'stabilityLLETest: SV >= 0');
pvt_assert(ischar(res.overall),        'LLE result.overall is a string');
pvt_assert(ischar(res.message),        'LLE result.message is a string');
pvt_assert(ischar(res.test1_message),  'LLE result.test1_message is a string');
pvt_assert(ischar(res.test2_message),  'LLE result.test2_message is a string');

pvt_assert(strcmp(res.overall, 'unstable'), ...
    'H2O/C10H22 300K/100bar: expected unstable (LLE)');
pvt_assert(any(flag == 2), ...
    'H2O/C10H22: at least one LLE test finds non-trivial solution');

% ---- Case 2: result.SL and result.SV match returned SL/SV ---------------
pvt_assert(res.SL == SL, 'LLE result.SL matches SL output');
pvt_assert(res.SV == SV, 'LLE result.SV matches SV output');

% ---- Case 3: FlashOptions accepted without error -------------------------
opts = FlashOptions();
opts.maxIteration = 80;
[flag3, ~, ~, res3] = stabilityLLETest(mix, thermo, opts);
pvt_assert(numel(flag3) == 2,    'stabilityLLETest with options: flag has 2 elements');
pvt_assert(ischar(res3.overall), 'stabilityLLETest with options: result.overall present');

% ---- Case 4: pure component → should be stable (trivial solution) --------
[comp_pure, ~] = addComponents({'C10H22'});
mix_pure = Mixture(comp_pure, 300, 100e5);
[flag4, ~, ~, res4] = stabilityLLETest(mix_pure, thermo);
pvt_assert(strcmp(res4.overall, 'stable'), ...
    'Pure C10H22: expected stable (single component)');
pvt_assert(all(flag4 == 1), ...
    'Pure C10H22: both LLE flags should be 1 (trivial)');

% ---- Case 5: dimethyl ether / water / decane (testcase3 system) ----------
[comp5, ~] = addComponents({'dimethyl ether', 'water', 'C10H22'});
mix5 = Mixture(comp5, 300, 100e5);
mix5.mole_fraction = [0.1 0.45 0.45];
mix5.bip.EOScons = [0 0.01 0.1; 0.01 0 0.1; 0.1 0.1 0];
[flag5, ~, ~, res5] = stabilityLLETest(mix5, thermo);
pvt_assert(numel(flag5) == 2, 'DME/H2O/C10H22: flag has 2 elements');
pvt_assert(any(strcmp(res5.overall, {'stable','unstable','inconclusive'})), ...
    'DME/H2O/C10H22: result.overall is valid');
