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
disp(P_bub)
disp(y_bub)
disp(flag_bp)