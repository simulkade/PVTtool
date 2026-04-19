function test_flash_vle()
% test_flash_vle  Integration tests for VLE flash.
%
% System: methanol (CH4O) + water (H2O) at 322.91 K with kij = -0.07.
% The negative flash traces the phase envelope using equimolar feed; vapor
% fractions outside [0,1] are expected when the feed is in the single-phase
% region (that is the purpose of the negative saturation method).
%
% We test:
%   1. vleflashnegative produces finite, normalised compositions.
%   2. At p = 26131 Pa (equimolar feed is two-phase) the vapor composition
%      matches experiment within 0.04 mole fraction.
%   3. Mass balance: z = V*y + (1-V)*x for a two-phase flash.
%   4. vleflash (standard) gives valid results in [0,1] for a two-phase feed.
%   5. kval_estimate produces positive K-values.

[component, ~] = addComponents({'CH4O', 'H2O'});
T0      = 322.91;          % [K]
thermo1 = addThermo();
thermo1.EOS = @PREOS;
options = FlashOptions();

% --- 1. Negative flash: compositions are finite and sum to 1 ---------------
% Three representative experimental pressures (Pa) from Ochi et al. (1986).
p_exp   = [15932, 26131, 52142];
y1_exp  = [0.2741, 0.6294, 0.9736];  % vapor methanol mole fraction

for i = 1:3
    mix = Mixture(component, T0, p_exp(i));
    mix.bip.EOScons = [0 -0.07; -0.07 0];
    [vapor_y, liquid_x, vapor_frac] = vleflashnegative(mix, thermo1, options);

    pvt_assert(isfinite(vapor_frac), ...
        sprintf('vleflashnegative p=%d: vapor_frac is finite', p_exp(i)));
    pvt_assert(isfinite(sum(vapor_y)) && isfinite(sum(liquid_x)), ...
        sprintf('vleflashnegative p=%d: compositions are finite', p_exp(i)));
    pvt_assert(abs(sum(vapor_y) - 1) < 1e-4, ...
        sprintf('vleflashnegative p=%d: vapor_y sums to 1', p_exp(i)));
    pvt_assert(abs(sum(liquid_x) - 1) < 1e-4, ...
        sprintf('vleflashnegative p=%d: liquid_x sums to 1', p_exp(i)));
end

% --- 2. Vapor composition matches experiment at middle pressure --------
% At p=26131 Pa with equimolar feed the vapor fraction is ~0.71 (two-phase);
% the vapor composition should agree with experiment within 0.04.
TOL_y = 0.04;
mix_mid = Mixture(component, T0, p_exp(2));
mix_mid.bip.EOScons = [0 -0.07; -0.07 0];
[vapor_y_mid, liquid_x_mid, vf_mid] = vleflashnegative(mix_mid, thermo1, options);

pvt_assert(vf_mid >= 0 && vf_mid <= 1, ...
    'vleflashnegative p=26131: equimolar feed is two-phase (vf in [0,1])');
pvt_assert(abs(vapor_y_mid(1) - y1_exp(2)) < TOL_y, ...
    sprintf('vleflashnegative p=26131: y(methanol) within %.2f of experiment', TOL_y));

% --- 3. Mass balance check (z = V*y + (1-V)*x) --------------------------
z = mix_mid.mole_fraction;
z_check = vf_mid * vapor_y_mid(1) + (1 - vf_mid) * liquid_x_mid(1);
pvt_assert(abs(z_check - z(1)) < 1e-3, 'vleflashnegative: mass balance satisfied');

% --- 4. Standard vleflash: two-phase feed with known compositions --------
% Use a feed composition (z1=0.35) that lies between x1=0.2131 and y1=0.6294
% at p=26131 Pa, guaranteeing a two-phase result.
mix_tp = Mixture(component, T0, 26131);
mix_tp.mole_fraction = [0.35 0.65];
mix_tp.bip.EOScons = [0 -0.07; -0.07 0];
[vy_tp, lx_tp, vf_tp] = vleflash(mix_tp, thermo1, options);

pvt_assert(vf_tp >= 0 && vf_tp <= 1,     'vleflash: vapor fraction in [0,1]');
pvt_assert(abs(sum(vy_tp) - 1) < 1e-6,   'vleflash: vapor_y sums to 1');
pvt_assert(abs(sum(lx_tp) - 1) < 1e-6,   'vleflash: liquid_x sums to 1');
pvt_assert(all(vy_tp >= 0),               'vleflash: vapor_y non-negative');
pvt_assert(all(lx_tp >= 0),               'vleflash: liquid_x non-negative');
% Check mass balance
z_tp = mix_tp.mole_fraction(1);
z_tp_check = vf_tp * vy_tp(1) + (1 - vf_tp) * lx_tp(1);
pvt_assert(abs(z_tp_check - z_tp) < 1e-4, 'vleflash: mass balance satisfied');

% --- 5. kval_estimate gives positive K-values ---------------------------
K = kval_estimate(mix_tp);
pvt_assert(numel(K) == 2, 'kval_estimate: returns N K-values');
pvt_assert(all(K > 0),    'kval_estimate: all K-values positive');
