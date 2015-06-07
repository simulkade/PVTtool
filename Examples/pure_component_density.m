% test case
% liquid-two-phase phase transition
clc; clear;
R=8.31415;
% Define the components and load pure physsical properties
[component, comp_flag] = addComponents({'H2O'});
% Define the thermodynamic models
T0 = 300; % [K]
p0 = 100e5; % [Pa]
thermo1 = addThermo();
thermo1.EOS = @PREOS;
mixture1 = addMixture(component, T0, p0);
% Define flash options
options.accuracy = 1e-7;
options.iteration = 100;
% call EOS to calculate z
T_range=linspace(300,750, 100);
vl=zeros(size(T_range));
vg=vl;
for i=1:length(T_range)
    mixture1.temperature=T_range(i);
    [liquid_z, vapor_z, fugacity, HR] = PREOS(mixture1, thermo1);
    vl(i)=liquid_z*R*T_range(i)/p0;
    vg(i)=vapor_z*R*T_range(i)/p0;
end
v_ig=R*T_range/p0;
steam_data=[293.15	1.79668542413313E-005
308.4393235246	1.8046698703882E-005
323.7286470491	1.81597323588409E-005
339.0179705737	1.83011794455744E-005
354.3072940982	1.8468517136899E-005
369.5966176228	1.86605872949479E-005
384.8859411473	1.88771933728929E-005
400.1752646719	1.91189243266625E-005
415.4645881964	1.93871042264311E-005
430.753911721	1.96838285872704E-005
446.0432352455	2.00120791010825E-005
461.3325587701	2.03759296011222E-005
476.6218822946	2.07808789544732E-005
491.9112058192	2.12343818951916E-005
507.2005293438	2.17467150701855E-005
522.4898528683	2.23324550333908E-005
537.7791763929	2.30131694106867E-005];
semilogx(vl, T_range, steam_data(:,2), steam_data(:,1), 'o')
xlabel('molar volume [m^3/mol]');
ylabel('Temperature [K]');
legend('Peng-Robinson', 'Steam Table');
