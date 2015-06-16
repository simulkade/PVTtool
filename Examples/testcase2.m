% test case
% liquid-two-phase phase transition
clc; clear;
% Define the components and load pure physsical properties
[component, comp_flag] = addComponents({'CH4', 'C2H6', 'C4H10', ...
'C6H14', 'C10H22', 'C15H32'});
% Define the thermodynamic models
T0 = 300; % [K]
p0 = 100e5; % [Pa]
thermo1 = addThermo();
thermo1.EOS = @PREOS;
mixture1 = addMixture(component, T0, p0);
% Define flash options
options.accuracy = 1e-7;
options.iteration = 100;
% Negative flash
[vapor_y, liquid_x, vapor_frac] = ...
        vleflashnegative(mixture1, thermo1, options)
mixture1.temperature=T_range(i);
[liquid_z, vapor_z, fugacity, HR] = PREOS(mixture1, thermo1);    

[s1, sl, sv] = stabilityTest(mixture1, thermo1)
