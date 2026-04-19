function test_residual_props()
% test_residual_props  Tests for residual property calculations.
%
% Tests SR, GR, VR, Cp_R, Cv_R from PREOS, PR78EOS, and SRKEOS.
% Physical sign checks and thermodynamic identity verification.
% System: pure CH4 (subcritical liquid at 150K/1MPa, vapor at 300K/5MPa).

[comp, ~] = addComponents({'CH4'});
R = 8.314;

% ---- PREOS subcritical liquid (CH4 150 K, 1 MPa) ----
mix_liq = Mixture(comp, 150, 1e6);
th = ThermoModel();
th.phase = 1; th.fugacity_switch = 0;
[zl, ~, ~, HR_l, props_l] = PREOS(mix_liq, th);

pvt_assert(isfinite(props_l.SR),   'PREOS liquid: SR is finite');
pvt_assert(isfinite(props_l.GR),   'PREOS liquid: GR is finite');
pvt_assert(isfinite(props_l.VR),   'PREOS liquid: VR is finite');
pvt_assert(isfinite(props_l.Cp_R), 'PREOS liquid: Cp_R is finite');
pvt_assert(isfinite(props_l.Cv_R), 'PREOS liquid: Cv_R is finite');

% Physical sign expectations for compressed liquid:
pvt_assert(props_l.SR < 0,   'PREOS liquid: SR < 0 (liquid more ordered than ideal gas)');
pvt_assert(props_l.VR < 0,   'PREOS liquid: VR < 0 (liquid denser than ideal gas)');
pvt_assert(props_l.HR < 0,   'PREOS liquid: HR < 0 (attractive interactions)');

% GR = HR - T*SR consistency
pvt_assert(abs(props_l.GR - (props_l.HR - 150*props_l.SR)) < 1e-6, ...
    'PREOS liquid: GR = HR - T*SR identity');

% VR = RT*(Z-1)/P
VR_check = R*150*(zl-1)/1e6;
pvt_assert(abs(props_l.VR - VR_check) < 1e-15, 'PREOS liquid: VR = RT(Z-1)/P');

% ---- PREOS vapor (CH4 300 K, 1 MPa — supercritical but gas-like) ----
mix_vap = Mixture(comp, 300, 1e6);
th.phase = 2;
[~, zv, ~, HR_v, props_v] = PREOS(mix_vap, th);

pvt_assert(props_v.VR > -1e-3,  'PREOS vapor: VR near zero or positive at low pressure');
VR_check_v = R*300*(zv-1)/1e6;
pvt_assert(abs(props_v.VR - VR_check_v) < 1e-15, 'PREOS vapor: VR = RT(Z-1)/P');

% ---- PR78EOS: same residual structure as PREOS for CH4 (omega=0.011 < 0.491) ----
th78 = ThermoModel(); th78.phase = 1; th78.fugacity_switch = 0;
[~, ~, ~, HR_78, props_78] = PR78EOS(mix_liq, th78);
pvt_assert(abs(HR_78 - HR_l)/abs(HR_l) < 0.01,      'PR78EOS: HR matches PREOS within 1% for CH4');
pvt_assert(abs(props_78.SR - props_l.SR)/abs(props_l.SR) < 0.01, ...
    'PR78EOS: SR matches PREOS within 1% for CH4');
pvt_assert(abs(props_78.GR - (props_78.HR - 150*props_78.SR)) < 1e-6, ...
    'PR78EOS: GR = HR - T*SR identity');

% ---- SRKEOS liquid (CH4 150K/1MPa) ----
th_srk = ThermoModel(); th_srk.EOS = @SRKEOS; th_srk.phase = 1; th_srk.fugacity_switch = 0;
[zl_srk, ~, ~, HR_srk, props_srk] = SRKEOS(mix_liq, th_srk);

pvt_assert(isfinite(props_srk.SR),   'SRKEOS: SR is finite');
pvt_assert(props_srk.SR < 0,         'SRKEOS liquid: SR < 0');
pvt_assert(props_srk.VR < 0,         'SRKEOS liquid: VR < 0');
pvt_assert(HR_srk < 0,               'SRKEOS liquid: HR < 0');

VR_srk_check = R*150*(zl_srk-1)/1e6;
pvt_assert(abs(props_srk.VR - VR_srk_check) < 1e-15, 'SRKEOS: VR = RT(Z-1)/P');
pvt_assert(abs(props_srk.GR - (props_srk.HR - 150*props_srk.SR)) < 1e-6, ...
    'SRKEOS: GR = HR - T*SR identity');

% ---- Cp_R > Cv_R for all EOS (always true for stable phases) ----
pvt_assert(props_l.Cp_R > props_l.Cv_R,   'PREOS: Cp_R > Cv_R (liquid)');
pvt_assert(props_srk.Cp_R > props_srk.Cv_R, 'SRKEOS: Cp_R > Cv_R (liquid)');

% ---- Two-component mixture: CH4 + C2H6, check finite props ----
[comp2, ~] = addComponents({'CH4', 'C2H6'});
mix2 = Mixture(comp2, 250, 5e6);
mix2.mole_fraction = [0.6 0.4];
th2 = ThermoModel(); th2.phase = 1; th2.fugacity_switch = 0;
[~, ~, ~, HR2, props2] = PREOS(mix2, th2);
pvt_assert(isfinite(HR2),          '2-comp PREOS: HR finite');
pvt_assert(isfinite(props2.SR),    '2-comp PREOS: SR finite');
pvt_assert(isfinite(props2.Cp_R),  '2-comp PREOS: Cp_R finite');
