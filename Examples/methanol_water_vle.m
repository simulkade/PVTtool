% test case
% liquid-two-phase phase transition
clc; clear;
% Define the components and load pure physsical properties
[component, comp_flag] = addComponents({'CH4O', 'H2O'});
% Define the thermodynamic models
T0 = 322.91; % [K]
p0 = [15.932
20.932
22.625
26.131
29.024
31.544
37.73
40.85
43.21
46.45
49.796
52.142]*1000; % [Pa]
x1=[0.0486
0.1218
0.1478
0.2131
0.2693
0.3252
0.5143
0.6279
0.7083
0.8037
0.9007
0.9461];
y1=[0.2741
0.4741
0.522
0.6294
0.7106
0.758
0.8203
0.8654
0.9007
0.9406
0.9627
0.9736];
thermo1 = addThermo();
thermo1.EOS = @PREOS;
x_calc=zeros(size(x1));
y_calc=x_calc;
mixture1 = addMixture(component, T0(1), p0(1));
BIP.EOScons = [0 -0.07; -0.07 0];
BIP.EOStdep = zeros(2);
mixture1.bip=BIP;
% Define flash options
options.accuracy = 1e-7;
options.iteration = 100;
% Negative flash
for i=1:length(x1)
  mixture1.pressure=p0(i);
  [vapor_y, liquid_x, vapor_frac] = ...
        vleflashnegative(mixture1, thermo1, options);
  x_calc(i)=liquid_x(1);
  y_calc(i)=vapor_y(1);
end
plot(p0, x1, 'o', p0, y1, 'o', p0, x_calc, p0, y_calc)  
 % [s1, sl, sv] = stabilityTest(mixture1, thermo1)   
