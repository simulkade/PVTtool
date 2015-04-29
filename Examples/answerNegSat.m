% test case
% liquid-two-phase phase transition
clc; clear;
% Define the components and load pure physsical properties
[component, comp_flag] = addComponents({'CH4', 'C2H6', 'C4H10', 'C6H14'});
% Define the thermodynamic models
thermo = struct();
thermo.mixingrule = 1; %van der Waals mixing rule
thermo.activity_model = @NRTL;
thermo.EOS = @PREOS;
thermo.phase = 1;
thermo.fugacity_switch = 1;
% define the mixture
BIP = struct();
BIP.EOScons = zeros(4);
BIP.EOStdep = zeros(4);
mixture = struct();
mixture.bip = BIP;
mixture.components = component;
mixture.mole_fraction = [0.3 0.25 0.35 0.1];
mixture.pressure = 100e5; %[Pa]
mixture.temperature = 300; % [K]

options = struct();
options.accuracy = 1e-7;
options.iteration = 100;

mixtureV = mixture;
mixtureL = mixture;
thermoV = thermo;
thermoL = thermo;
thermoL.phase = 1;
thermoV.phase = 2;
nData = 100;
p = linspace(60,80, nData);
T = 300; % [K]
R = 8.314; % J/(mol.k)
liquid_x = zeros(nData, 4);
vapor_y = zeros(nData, 4);
vapor_frac = zeros(nData, 1);
Cg = zeros(nData, 1);
Cl = zeros(nData, 1);
for i = 1:nData
    mixture.pressure = p(i)*1e5; %[Pa]
    [vapor_y(i,:), liquid_x(i,:), vapor_frac(i)] = ...
        vleflashnegative(mixture, thermo, options);
    mixtureV.mole_fraction = vapor_y(i,:);
    mixtureV.pressure = p(i)*1e5;
    mixtureL.mole_fraction = liquid_x(i,:);
    mixtureL.pressure = p(i)*1e5;
    liqz = PREOS(mixtureL,thermoL);
    vapz = PREOS(mixtureV,thermoV);
    Cg(i) = p(i)*1e5/(R*T*vapz);
    Cl(i) = p(i)*1e5/(R*T*liqz); 
end
Cg(vapor_frac<=0) = Cl(vapor_frac<=0);
Sg = (vapor_frac./Cg)./(vapor_frac./Cg+(1-vapor_frac)./Cl);

%% visualize
figure(1);
subplot(2,1,1);
plot(p, Sg);
subplot(2,1,2);
plot(p, (Sg>0).*Sg);
xlabel('pressure [bar]');
ylabel('Saturation');


figure(2);
subplot(3,1,1);
plot(p, liquid_x);
ylabel('liquid mole fraction');
subplot(3,1,2);
plot(p, vapor_y);
ylabel('vapor mole fraction');
subplot(3,1,3);
vf_rep = repmat(vapor_frac, 1, 4);
xl = (vf_rep<=0).*(vf_rep.*vapor_y+(1-vf_rep).*liquid_x)+(vf_rep>0).*liquid_x;
plot(p, xl)
ylabel('vapor mole fraction');
xlabel('pressure [bar]');
legend('CH4', 'C2H6', 'C4H10', 'C6H14');