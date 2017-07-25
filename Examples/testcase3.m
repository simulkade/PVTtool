% test case
% liquid-two-phase phase transition
clc; clear;
% Define the components and load pure physsical properties
[component, comp_flag] = addComponents({'dimethyl ether', 'water', 'C10H22'});
% Define the thermodynamic models
T0 = 300; % [K]
p0 = 100e5; % [Pa]
thermo1 = addThermo();
thermo1.EOS = @PREOS;
mixture1 = addMixture(component, T0, p0);
mixture1.mole_fraction = [0.1, 0.45, 0.45];
BIP = struct();
BIP.EOScons = [0 0.01 0.1; 0.01 0 0.1; 0.1 0.1 0];
BIP.EOStdep = zeros(3);
mixture1.bip = BIP;
% Define flash options
options.accuracy = 1e-7;
options.iteration = 100;
% Negative flash
[vapor_y, liquid_x, vapor_frac] = ...
        lleflash(mixture1, thermo1, options)
%mixture1.temperature=T0;
[liquid_z, vapor_z, fugacity, HR] = PREOS(mixture1, thermo1);    

[s1, sl, sv] = stabilityLLETest(mixture1, thermo1)
