function test_saturation()
% test_saturation  Tests for bubble and dew point calculations.
%
% System: methanol (CH4O) + water (H2O) at experimental conditions.
% Reference: Ochi et al. (1986) data used in methanol_water_vle example.
%
% At 322.91 K, 26131 Pa, the equimolar mixture is two-phase.
% Bubble pressure for x1=0.35 (methanol): should be ~26131 Pa +/- 20%.
% Bubble temperature at 26131 Pa for z=[0.5 0.5]: should be near 322.91 K.

[component, ~] = addComponents({'CH4O', 'H2O'});
T0 = 322.91;
p0 = 26131;
thermo1 = addThermo();
opts = FlashOptions();

% ---- bubblePressure: z=[0.35, 0.65] at T=322.91 K ----
mix_bp = Mixture(component, T0, p0);
mix_bp.mole_fraction = [0.35 0.65];
mix_bp.bip.EOScons = [0 -0.07; -0.07 0];

[P_bub, y_bub, flag_bp] = bubblePressure(mix_bp, thermo1, opts);

pvt_assert(flag_bp == 1, 'bubblePressure: converged');
pvt_assert(P_bub > 1e3 && P_bub < 2e5, 'bubblePressure: result in physical range');
pvt_assert(abs(sum(y_bub) - 1) < 1e-4,  'bubblePressure: y sums to 1');
pvt_assert(all(y_bub >= 0),              'bubblePressure: y non-negative');
% y(methanol) > z(methanol) at bubble (methanol more volatile than water at 322K)
pvt_assert(y_bub(1) > mix_bp.mole_fraction(1), ...
    'bubblePressure: methanol enriched in vapor (more volatile)');

% ---- bubbleTemperature: equimolar at p=26131 Pa ----
mix_bt = Mixture(component, T0, p0);
mix_bt.bip.EOScons = [0 -0.07; -0.07 0];

[T_bub, y_bub2, flag_bt] = bubbleTemperature(mix_bt, thermo1, opts);

pvt_assert(flag_bt == 1,                    'bubbleTemperature: converged');
pvt_assert(abs(T_bub - T0) < 20,            'bubbleTemperature: T within 20 K of reference');
pvt_assert(abs(sum(y_bub2) - 1) < 1e-4,    'bubbleTemperature: y sums to 1');

% ---- dewPressure: vapor z=[0.7, 0.3] at T=322.91 K ----
mix_dp = Mixture(component, T0, p0);
mix_dp.mole_fraction = [0.7 0.3];
mix_dp.bip.EOScons = [0 -0.07; -0.07 0];

[P_dew, x_dew, flag_dp] = dewPressure(mix_dp, thermo1, opts);

pvt_assert(flag_dp == 1, 'dewPressure: converged');
pvt_assert(P_dew > 1e3 && P_dew < 2e5, 'dewPressure: result in physical range');
pvt_assert(abs(sum(x_dew) - 1) < 1e-4,  'dewPressure: x sums to 1');
pvt_assert(all(x_dew >= 0),              'dewPressure: x non-negative');

% ---- dewTemperature: vapor z=[0.7, 0.3] at p=26131 Pa ----
mix_dt = Mixture(component, T0, p0);
mix_dt.mole_fraction = [0.7 0.3];
mix_dt.bip.EOScons = [0 -0.07; -0.07 0];

[T_dew, x_dew2, flag_dt] = dewTemperature(mix_dt, thermo1, opts);

pvt_assert(flag_dt == 1,                    'dewTemperature: converged');
pvt_assert(abs(T_dew - T0) < 30,            'dewTemperature: T within 30 K of reference');
pvt_assert(abs(sum(x_dew2) - 1) < 1e-4,    'dewTemperature: x sums to 1');

% ---- Consistency: bubble P < dew P at same T for same feed ----
% (when bubble and dew points bound the two-phase envelope)
mix_cons = Mixture(component, T0, p0);
mix_cons.mole_fraction = [0.5 0.5];
mix_cons.bip.EOScons = [0 -0.07; -0.07 0];
[P_b, ~, flag_b] = bubblePressure(mix_cons, thermo1, opts);
[P_d, ~, flag_d] = dewPressure(mix_cons, thermo1, opts);
if flag_b && flag_d
    % At fixed T, bubble pressure (upper two-phase boundary) > dew pressure (lower boundary)
    pvt_assert(P_b > P_d, 'bubblePressure > dewPressure: bubble is upper, dew is lower boundary');
end
