% test case
% liquid-two-phase phase transition
clc; clear;
% Define the components and load pure physsical properties
z=[0.1, 0.1, 0.1, 0.2, 0.3, 0.2]
[component, comp_flag] = addComponents({'CH4', 'C2H6', 'C4H10', ...
'C6H14', 'C8H18', 'C10H22'});
% Define the thermodynamic models
T = 100+273.15; % [K]
T0= 25+273.15; % [K] surface temperature
p0= 1e5; % [Pa] surface pressure
thermo1 = addThermo();
thermo1.EOS = @PREOS;
mixture1 = addMixture(component, T0, p0);
mixture1.mole_fraction = z;
% Define flash options
options.accuracy = 1e-7;
options.iteration = 100;
% flash for surface properties
[vapor_y, liquid_x, vapor_frac_surf] = ...
        vleflash(mixture1, thermo1, options)
[liquid_z, vapor_z, fugacity, HR] = PREOS(mixture1, thermo1);
R=8.314; % [J/(mol.K)]
mixture1.mole_fraction=liquid_x;
[liquid_z_surf, vapor_z, fugacity, HR] = PREOS(mixture1, thermo1);
V_oil_surface = liquid_z_surf*(1-vapor_frac_surf)*R*T0/p0;
mixture1.mole_fraction=vapor_y;
[liquid_z, vapor_z_surf, fugacity, HR] = PREOS(mixture1, thermo1);
V_gas_surface = vapor_z_surf*vapor_frac_surf*R*T0/p0;

n_point=50;
p_range=linspace(p0, 300e5, n_point);
for i=1:n_point
p=p_range(i);
% flash for reservoir properties
mixture1.temperature=T;
mixture1.pressure=p;
mixture1.mole_fraction = z;
[vapor_y, liquid_x, vapor_frac_res] = ...
        vleflashnegative(mixture1, thermo1, options)
if vapor_frac_res<=0
  vapor_y=z;
  liquid_x=z;
  vapor_frac_res=0;
end  
mixture1.mole_fraction=liquid_x;
[liquid_z_res, vapor_z, fugacity, HR] = PREOS(mixture1, thermo1);
V_oil_res = liquid_z_res*(1-vapor_frac_res)*R*T/p;
mixture1.mole_fraction=vapor_y;
[liquid_z, vapor_z_res, fugacity, HR] = PREOS(mixture1, thermo1);
V_gas_res = vapor_z_res*vapor_frac_res*R*T/p;

%[s1, sl, sv] = stabilityTest(mixture1, thermo1)
Bo(i) = V_oil_res/V_oil_surface
Bg(i) = V_gas_res/V_oil_surface
Rs(i) = (vapor_frac_surf-vapor_frac_res)*vapor_z_surf*R*T0/p0/V_oil_surface
rs1 = (vapor_frac_res-vapor_frac_surf)*liquid_z_surf*R*T0/p0/V_gas_surface;
if rs1<0
  rs1=0;
end
rs(i)=rs1;
end