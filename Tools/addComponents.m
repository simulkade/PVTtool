function [component, flag] = addComponents(component_name_formula)
% This function returns an structure array of components including the
% physical properties of pure components. The argument
% component_name_formula is defined as:
% component_name_formula = {'comp1_name_or_formula',
% 'comp2_name_or_formula', ...}
% if one of the components does not exist, the function returns the index to that component in flag
%
% SYNOPSIS:
%
%
% PARAMETERS:
%
%
% RETURNS:
%
%
% EXAMPLE:
%
% SEE ALSO:
%

%{
Copyright (c) 2012, 2013, Ali Akbar Eftekhari
All rights reserved.

Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following
conditions are met:

    *   Redistributions of source code must retain the above copyright notice,
        this list of conditions and the following disclaimer.
    *   Redistributions in binary form must reproduce the above
        copyright notice, this list of conditions and the following
        disclaimer in the documentation and/or other materials provided
        with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%}

component = struct([]); % define an array of component structure
N = length(component_name_formula);
data_index = zeros(1,N);
for n = 1:N
     temp_index = find(strcmpi(component_formula,component_name_formula(n)),1);
     if (isempty(temp_index)==0)
         data_index(n) = temp_index;
     else
         temp_index = find(strcmpi(component_name,component_name_formula(n)),1);
         if (isempty(temp_index)==0)
             data_index(n) = temp_index;
         else
            data_index(n) = 0;
         end
     end
end
flag = find(data_index == 0);
data_index = data_index(data_index>0);
N = length(data_index);
for n = 1:N
    comp_index = data_index(n);
    component(n).name = component_name(comp_index);
    component(n).formula = component_formula(comp_index);
    component(n).MW = Molecular_weight_data(comp_index);
    component(n).Tc = critical_temperature_data(comp_index);
    component(n).Pc = critical_pressure_data(comp_index);
    component(n).Vc = critical_volume_data(comp_index);
    component(n).Zc = critical_compressibility_factor(comp_index);
    component(n).acentric_factor = acentric_factor(comp_index);
    component(n).Psat_eq = vapor_pressure_equation;
    component(n).Psat_coefs = vapor_pressure_coefs(comp_index,:);
    component(n).PsatTrange = vapor_pressure_T_range(comp_index,:);
    component(n).Psatrange = vapor_pressure_p_range(comp_index,:);
    component(n).dh_vap_eq = dh_vaporization_equation;
    component(n).dh_vap_coefs = dh_vaporization_coefs(comp_index,:);
    component(n).dh_vap_Trange = dh_vaporization_T_range(comp_index,:);
    component(n).dh_vap_range = dh_vaporization_dhv_range(comp_index,:);
    if liquid_heat_capacity_equation_num(comp_index)==1
        component(n).cp_liq_eq = liquid_heat_capacity_equation_1;
    else
        component(n).cp_liq_eq = liquid_heat_capacity_equation_2;
    end
    component(n).cp_liq_coefs = liquid_heat_capacity_coef(comp_index,:);
    component(n).cp_liq_Trange = liquid_heat_capacity_T_range(comp_index,:);
    component(n).cp_liq_range = liquid_heat_capacity_cp_range(comp_index,:);
    if ideal_gas_heat_capacity_equation_num(comp_index)==1
        component(n).cp_ig_eq = ideal_gas_heat_capacity_equation_1;
    elseif ideal_gas_heat_capacity_equation_num(comp_index)==2
        component(n).cp_ig_eq = ideal_gas_heat_capacity_equation_2;
    else
        component(n).cp_ig_eq = ideal_gas_heat_capacity_equation_3;
    end
    component(n).cp_ig_coefs = ideal_gas_heat_capacity_coefs(comp_index,:);
    component(n).cp_ig_Trange = ideal_gas_heat_capacity_T_range(comp_index,:);
    component(n).cp_ig_range = ideal_gas_heat_capacity_cp_range(comp_index,:);
    component(n).dhf_ig = dh_formation_ig(comp_index);
    component(n).dgf_ig = dg_formation_ig(comp_index);
    component(n).ds_ig = ds_ideal_gas(comp_index);
    component(n).dh_comb = dh_combustion(comp_index);
end
