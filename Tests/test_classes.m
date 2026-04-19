function test_classes()
% test_classes  Unit tests for BIP, Component, Mixture, ThermoModel, FlashOptions.

% --- BIP ---
b2 = BIP(2);
pvt_assert(isequal(size(b2.EOScons),  [2 2]), 'BIP(2): EOScons is 2x2');
pvt_assert(isequal(size(b2.NRTLalfa), [2 2]), 'BIP(2): NRTLalfa is 2x2');
pvt_assert(isequal(size(b2.UNIQUACR), [1 2]), 'BIP(2): UNIQUACR is 1x2');
pvt_assert(all(b2.EOScons(:) == 0),           'BIP(2): EOScons initialised to 0');
pvt_assert(all(b2.NRTLcons(:) == 0),          'BIP(2): NRTLcons initialised to 0');

b3 = BIP(3);
pvt_assert(isequal(size(b3.UNIQUACcons), [3 3]), 'BIP(3): UNIQUACcons is 3x3');

% --- Component.fromDatabaseArray ---
[comps, flag] = Component.fromDatabaseArray({'CH4', 'C2H6'});
pvt_assert(isempty(flag),                    'addComponents: no missing components');
pvt_assert(numel(comps) == 2,               'addComponents: returns 2 components');
pvt_assert(isa(comps(1), 'Component'),      'addComponents: returns Component objects');
pvt_assert(~isempty(comps(1).name),         'Component: name is populated');
pvt_assert(comps(1).Tc > 0,                 'Component: Tc > 0');
pvt_assert(comps(1).Pc > 0,                 'Component: Pc > 0');
pvt_assert(comps(1).acentric_factor >= 0,   'Component: acentric factor >= 0');

% --- Component.fromDatabase (single) ---
c = Component.fromDatabase('H2O');
pvt_assert(isa(c, 'Component'), 'fromDatabase: returns Component');
pvt_assert(c.Tc > 0,            'fromDatabase: Tc populated');

% --- addComponents wrapper ---
[comps2, flag2] = addComponents({'CH4', 'H2O'});
pvt_assert(isempty(flag2),         'addComponents wrapper: no missing');
pvt_assert(numel(comps2) == 2,    'addComponents wrapper: 2 components');
pvt_assert(isa(comps2(1), 'Component'), 'addComponents wrapper: Component type');

% --- Missing component ---
[~, flag_bad] = addComponents({'ThisDoesNotExist999'});
pvt_assert(~isempty(flag_bad), 'addComponents: missing component flag is non-empty');

% --- Mixture ---
[comp, ~] = addComponents({'CH4', 'C2H6'});
mix = Mixture(comp, 300, 1e6);
pvt_assert(isa(mix, 'Mixture'),           'Mixture: correct type');
pvt_assert(mix.temperature == 300,        'Mixture: temperature set');
pvt_assert(mix.pressure == 1e6,           'Mixture: pressure set');
pvt_assert(numel(mix.mole_fraction) == 2, 'Mixture: mole_fraction length');
pvt_assert(abs(sum(mix.mole_fraction)-1) < 1e-12, 'Mixture: mole fractions sum to 1');
pvt_assert(isa(mix.bip, 'BIP'),           'Mixture: bip is BIP object');

% --- addMixture wrapper ---
mix2 = addMixture(comp, 350, 2e6);
pvt_assert(isa(mix2, 'Mixture'), 'addMixture: returns Mixture');
pvt_assert(mix2.temperature == 350, 'addMixture: temperature set');

% --- ThermoModel ---
th = ThermoModel();
pvt_assert(isa(th, 'ThermoModel'),     'ThermoModel: correct type');
pvt_assert(th.mixingrule == 1,         'ThermoModel: default mixing rule = 1');
pvt_assert(th.fugacity_switch == 1,    'ThermoModel: default fugacity switch = 1');
pvt_assert(isequal(th.EOS, @PREOS),   'ThermoModel: default EOS is PREOS');

% --- addThermo wrapper ---
th2 = addThermo();
pvt_assert(isa(th2, 'ThermoModel'), 'addThermo: returns ThermoModel');

% --- FlashOptions ---
opts = FlashOptions();
pvt_assert(isa(opts, 'FlashOptions'),     'FlashOptions: correct type');
pvt_assert(opts.accuracy == 1e-7,         'FlashOptions: default accuracy');
pvt_assert(opts.iteration == 100,         'FlashOptions: default iteration');
pvt_assert(opts.convergenceMaxError==1e-10,'FlashOptions: default convergenceMaxError');
