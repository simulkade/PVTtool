% test case
clc; clear;
% Define the components and load pure physsical properties
[component, comp_flag] = addComponents({'CO2', 'C10H22'});
% Define the thermodynamic models
thermo = struct();
thermo.mixingrule = 1; %van der Waals mixing rule
thermo.activity_model = @NRTL;
thermo.EOS = @PREOS;
thermo.phase = 1;
thermo.fugacity_switch = 1;
% define the mixture
% -------------------------------
% BIP data CO2-n-Decane
% Tdata = [462 477 543 584];
% bipco2decane = [0.1058 0.1182 0.1735 0.2856];
% bipdata = polyfit(Tdata(1:end-1), bipco2decane(1:end-1), 1);
T_exp = [277.6 310.9 344.3 377.6 410.9 444.3 477.6 510.9];
bip_pr = [0.1221 0.1106 0.1039 0.1020 0.1049 0.1127 0.1254 0.1428];
bipT = @(T)(interp1(T_exp, bip_pr, T));
% --------------------------------
R = 8.314472; % J/mol/K
T = 312; % K
i = 0;
M = [component.MW]; % kg/mol
for p = 10e5:10e5:70e5; % Pa
    i = i+1;
    BIP = struct();
    BIP.EOScons = [0 bipT(T); bipT(T) 0];
    BIP.EOStdep = zeros(2);
    mixture = struct();
    mixture.bip = BIP;
    mixture.components = component;
    mixture.mole_fraction = [0.9 0.1];
    mixture.pressure = p; %[Pa]
    mixture.temperature = T; % [K]

    options = struct();
    options.accuracy = 1e-7;
    options.iteration = 200;
    [vapor_y, liquid_x, vapor_frac(i)] = vleflash(mixture, thermo, options);
    x_CO2(i) = liquid_x(1);
    mixture2 = mixture;
    mixture2.mole_fraction = liquid_x;
    liqz = PREOS(mixture2,thermo);
    c(i)= p/(R*T*liqz); % mol/m^3
    H(i) = p/(c(i)*liquid_x(1));
    MW = sum(liquid_x.*M);
    rho(i) = c(i)*MW;
end

mixture.pressure = p; %[Pa]
[vapor_y, liquid_x, vapor_frac(i)] = vleflash(mixture, thermo, options);
i=0;
for x = 0:(liquid_x(1)/10):liquid_x(1)
    i=i+1;
    mixture2.mole_fraction = [x 1-x];
    liqz = PREOS(mixture2,thermo);
    cc(i)= p/(R*T*liqz); % mol/m^3
    MW = sum([x 1-x].*M);
    rhom(i) = cc(i)*MW;
    
end
xco2 = 0:(liquid_x(1)/10):liquid_x(1);
pi = 10e5:10e5:70e5;

A = [pi' H'];
B = [(cc.*xco2)' rhom'];
%xlswrite('henry.xlsx', A);
%xlswrite('density.xlsx', B);
plot(pi, H, 'o', pi,H);
xlabel('Pressure (Pa)');
ylabel('Henry''s constant [Pa.m^3/mol]'); 
plot(pi, x_CO2, 'o', pi,x_CO2);
plot(pi, c, 'o', pi,c);
% xlswrite('co2decane', [pi'