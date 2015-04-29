function vf = testflashpT(p, T)
% test case
% liquid-two-phase phase transition
% clc; clear;
% Define the components and load pure physsical properties
[component, ~] = addComponents({'CH4', 'C2H6', 'C4H10', 'C6H14'});
% Define the thermodynamic models
T0 = T; % [K]
p0 = p; % [Pa]
thermo1 = addThermo();
mixture1 = addMixture(component, T0, p0);
% Define flash options
options.accuracy = 1e-7;
options.iteration = 100;
% Negative flash
[~, ~, vf] = ...
        vleflashnegative(mixture1, thermo1, options);
