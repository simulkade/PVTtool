function test_eos()
% test_eos  Verify PREOS, SRKEOS, and PR78EOS Z-factors and fugacities.
%
% Reference conditions: pure CH4 at 150 K, 1 MPa (subcritical: Tc=190.6 K,
% Pc=4.6 MPa) — guarantees two distinct real roots (liquid and vapor Z).
% Supercritical tests (300 K, 10 MPa) use single-root assertions only.

[comp, ~] = addComponents({'CH4'});

% ---- Subcritical conditions: two distinct roots expected ----
T_sub = 150;    % [K]  below Tc = 190.6 K
p_sub = 1e6;    % [Pa]  below Pc = 4.6 MPa
mix_sub = Mixture(comp, T_sub, p_sub);
th_sub = ThermoModel();
th_sub.fugacity_switch = 1;

th_sub.phase = 1;
[zl, zv, fug_l, HR] = PREOS(mix_sub, th_sub);
pvt_assert(zl > 0,    'PREOS subcritical: Zl > 0');
pvt_assert(zv > 0,    'PREOS subcritical: Zv > 0');
pvt_assert(zv > zl,   'PREOS subcritical: Zv > Zl (two-phase region)');
pvt_assert(zl < 0.1,  'PREOS subcritical: Zl in liquid range');
pvt_assert(zv > 0.7,  'PREOS subcritical: Zv in vapor range');
pvt_assert(fug_l > 0, 'PREOS subcritical: liquid fugacity > 0');

th_sub.phase = 2;
[~, ~, fug_v, ~] = PREOS(mix_sub, th_sub);
pvt_assert(fug_v > 0, 'PREOS subcritical: vapor fugacity > 0');

% K-value from liquid/vapor fugacity ratio
Kval = fug_l / fug_v;
pvt_assert(Kval > 0 && Kval < 1e6, 'PREOS: K-value in reasonable range');

% ---- Supercritical: only one real root, Zl == Zv ----
T_sup = 300;    % [K]  above Tc = 190.6 K
p_sup = 10e6;   % [Pa]  above Pc = 4.6 MPa → supercritical
mix_sup = Mixture(comp, T_sup, p_sup);
th_sup = ThermoModel();
th_sup.phase = 1; th_sup.fugacity_switch = 1;
[zl_sc, zv_sc, fug_sc, ~] = PREOS(mix_sup, th_sup);
pvt_assert(zl_sc > 0,         'PREOS supercritical: Z > 0');
pvt_assert(zl_sc == zv_sc,    'PREOS supercritical: Zl == Zv (single root)');
pvt_assert(fug_sc > 0,        'PREOS supercritical: fugacity > 0');

% ---- SRKEOS (subcritical) ----
[zl_srk, zv_srk, fug_srk, ~] = SRKEOS(mix_sub, th_sub);
pvt_assert(zl_srk > 0,  'SRKEOS: Zl > 0');
pvt_assert(zv_srk > 0,  'SRKEOS: Zv > 0');
pvt_assert(fug_srk > 0, 'SRKEOS: fugacity > 0');

% ---- PR78EOS matches PREOS for CH4 (acentric factor 0.011 < 0.491) ----
th_sub.phase = 1;
[zl_pr78, zv_pr78, ~, ~] = PR78EOS(mix_sub, th_sub);
pvt_assert(abs(zl_pr78 - zl) / zl < 0.01, 'PR78EOS: Zl matches PREOS within 1% for CH4');
pvt_assert(abs(zv_pr78 - zv) / zv < 0.01, 'PR78EOS: Zv matches PREOS within 1% for CH4');

% ---- Pure water (high acentric factor 0.345 > 0.491, so PR78 differs) ----
[comp_w, ~] = addComponents({'H2O'});
mix_w = Mixture(comp_w, 400, 1e6);
th2 = ThermoModel();
th2.phase = 2; th2.fugacity_switch = 1;
[~, zv_w, fug_w, ~] = PREOS(mix_w, th2);
pvt_assert(zv_w > 0.9, 'PREOS water vapor: Zv close to ideal (low pressure)');
pvt_assert(fug_w > 0,  'PREOS water vapor: fugacity > 0');

% ---- Residual enthalpy has correct sign (negative for liquid) ----
th3 = ThermoModel();
th3.phase = 1; th3.fugacity_switch = 0;
[~, ~, ~, HR_liq] = PREOS(mix_sub, th3);
pvt_assert(HR_liq < 0, 'PREOS: residual enthalpy negative for liquid CH4');

% ---- Two-component mixture: PREOS with van der Waals mixing ----
[comp2, ~] = addComponents({'CH4', 'C2H6'});
mix2 = Mixture(comp2, 250, 5e6);
mix2.mole_fraction = [0.6 0.4];
th4 = ThermoModel();  % vdW, phase=1
[zl2, zv2, fug2, ~] = PREOS(mix2, th4);
pvt_assert(zl2 > 0 && zv2 > 0, 'PREOS 2-comp: both Z roots positive');
pvt_assert(numel(fug2) == 2,   'PREOS 2-comp: fugacity has 2 elements');
pvt_assert(all(fug2 > 0),       'PREOS 2-comp: all fugacities positive');
