classdef Component
% Component  Pure-component thermodynamic properties.
%
%   Typical usage: load from database
%     comp = Component.fromDatabase('CH4');
%     comps = Component.fromDatabaseArray({'CH4', 'C2H6'});
%
%   Or construct directly:
%     comp = Component(name, formula, MW, Tc, Pc, Vc, Zc,
%                      acentric_factor, Psat_eq, Psat_coefs,
%                      PsatTrange, Psatrange,
%                      dh_vap_eq, dh_vap_coefs, dh_vap_Trange, dh_vap_range,
%                      cp_liq_eq, cp_liq_coefs, cp_liq_Trange, cp_liq_range,
%                      cp_ig_eq,  cp_ig_coefs,  cp_ig_Trange,  cp_ig_range,
%                      dhf_ig, dgf_ig, ds_ig, dh_comb)

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

    properties
        name
        formula
        MW
        Tc
        Pc
        Vc
        Zc
        acentric_factor
        Psat_eq
        Psat_coefs
        PsatTrange
        Psatrange
        dh_vap_eq
        dh_vap_coefs
        dh_vap_Trange
        dh_vap_range
        cp_liq_eq
        cp_liq_coefs
        cp_liq_Trange
        cp_liq_range
        cp_ig_eq
        cp_ig_coefs
        cp_ig_Trange
        cp_ig_range
        dhf_ig
        dgf_ig
        ds_ig
        dh_comb
    end

    methods
        function obj = Component(name, formula, ...
                MW, Tc, Pc, Vc, Zc, acentric_factor, ...
                Psat_eq, Psat_coefs, PsatTrange, Psatrange, ...
                dh_vap_eq, dh_vap_coefs, dh_vap_Trange, dh_vap_range, ...
                cp_liq_eq, cp_liq_coefs, cp_liq_Trange, cp_liq_range, ...
                cp_ig_eq, cp_ig_coefs, cp_ig_Trange, cp_ig_range, ...
                dhf_ig, dgf_ig, ds_ig, dh_comb)
            % Zero-argument call is used by Octave/MATLAB to pre-allocate
            % object arrays (e.g. obj(n) = Component(...)).
            if nargin == 0; return; end
            obj.name            = name;
            obj.formula         = formula;
            obj.MW              = MW;
            obj.Tc              = Tc;
            obj.Pc              = Pc;
            obj.Vc              = Vc;
            obj.Zc              = Zc;
            obj.acentric_factor = acentric_factor;
            obj.Psat_eq         = Psat_eq;
            obj.Psat_coefs      = Psat_coefs;
            obj.PsatTrange      = PsatTrange;
            obj.Psatrange       = Psatrange;
            obj.dh_vap_eq       = dh_vap_eq;
            obj.dh_vap_coefs    = dh_vap_coefs;
            obj.dh_vap_Trange   = dh_vap_Trange;
            obj.dh_vap_range    = dh_vap_range;
            obj.cp_liq_eq       = cp_liq_eq;
            obj.cp_liq_coefs    = cp_liq_coefs;
            obj.cp_liq_Trange   = cp_liq_Trange;
            obj.cp_liq_range    = cp_liq_range;
            obj.cp_ig_eq        = cp_ig_eq;
            obj.cp_ig_coefs     = cp_ig_coefs;
            obj.cp_ig_Trange    = cp_ig_Trange;
            obj.cp_ig_range     = cp_ig_range;
            obj.dhf_ig          = dhf_ig;
            obj.dgf_ig          = dgf_ig;
            obj.ds_ig           = ds_ig;
            obj.dh_comb         = dh_comb;
        end

        function p_sat = vapor_pressure(obj, T)
            % vapor_pressure(T)  returns p_sat [Pa] at temperature T [K]
            p_sat = obj.Psat_eq(T, obj.Psat_coefs);
        end

        function dh_v = dh_vap(obj, T)
            % dh_vap(T)  returns enthalpy of vaporisation [J/mol] at T [K]
            dh_v = obj.dh_vap_eq(T/obj.Tc, obj.dh_vap_coefs);
        end

        function cp_l = cp_liq(obj, T)
            % cp_liq(T)  returns liquid heat capacity [J/(mol·K)] at T [K]
            cp_l = obj.cp_liq_eq(T/obj.Tc, obj.cp_liq_coefs);
        end

        function cp = cp_ig(obj, T)
            % cp_ig(T)  returns ideal-gas heat capacity [J/(mol·K)] at T [K]
            cp = obj.cp_ig_eq(T, obj.cp_ig_coefs);
        end
    end

    methods (Static)
        function comp = fromDatabase(name_or_formula)
            % fromDatabase  Load a single Component from puredata.mat.
            %   comp = Component.fromDatabase('CH4')
            %   comp = Component.fromDatabase('methane')
            comps = Component.fromDatabaseArray({name_or_formula});
            if isempty(comps)
                error('Component:notFound', ...
                    'Component "%s" not found in database.', name_or_formula);
            end
            comp = comps(1);
        end

        function [components, flag] = fromDatabaseArray(name_formula_cell)
            % fromDatabaseArray  Load multiple Components from puredata.mat.
            %   [comps, flag] = Component.fromDatabaseArray({'CH4','C2H6'})
            %   flag contains indices of names that were not found.
            load puredata.mat
            N = length(name_formula_cell);
            data_index = zeros(1, N);
            for n = 1:N
                idx = find(strcmpi(component_formula, name_formula_cell(n)), 1);
                if ~isempty(idx)
                    data_index(n) = idx;
                else
                    idx = find(strcmpi(component_name, name_formula_cell(n)), 1);
                    if ~isempty(idx)
                        data_index(n) = idx;
                    end
                end
            end
            flag = find(data_index == 0);
            data_index = data_index(data_index > 0);
            Nfound = length(data_index);
            if Nfound == 0
                components = [];
                return
            end
            % Pre-allocate by constructing the last element first so MATLAB
            % knows the array type, then fill in order.
            i = Nfound;
            ci = data_index(i);
            if liquid_heat_capacity_equation_num(ci) == 1
                cp_liq_fn = liquid_heat_capacity_equation_1;
            else
                cp_liq_fn = liquid_heat_capacity_equation_2;
            end
            if ideal_gas_heat_capacity_equation_num(ci) == 1
                cp_ig_fn = ideal_gas_heat_capacity_equation_1;
            elseif ideal_gas_heat_capacity_equation_num(ci) == 2
                cp_ig_fn = ideal_gas_heat_capacity_equation_2;
            else
                cp_ig_fn = ideal_gas_heat_capacity_equation_3;
            end
            components(Nfound) = Component( ...
                component_name{ci}, component_formula{ci}, ...
                Molecular_weight_data(ci), ...
                critical_temperature_data(ci), critical_pressure_data(ci), ...
                critical_volume_data(ci), critical_compressibility_factor(ci), ...
                acentric_factor(ci), ...
                vapor_pressure_equation, vapor_pressure_coefs(ci,:), ...
                vapor_pressure_T_range(ci,:), vapor_pressure_p_range(ci,:), ...
                dh_vaporization_equation, dh_vaporization_coefs(ci,:), ...
                dh_vaporization_T_range(ci,:), dh_vaporization_dhv_range(ci,:), ...
                cp_liq_fn, liquid_heat_capacity_coef(ci,:), ...
                liquid_heat_capacity_T_range(ci,:), liquid_heat_capacity_cp_range(ci,:), ...
                cp_ig_fn, ideal_gas_heat_capacity_coefs(ci,:), ...
                ideal_gas_heat_capacity_T_range(ci,:), ideal_gas_heat_capacity_cp_range(ci,:), ...
                dh_formation_ig(ci), dg_formation_ig(ci), ...
                ds_ideal_gas(ci), dh_combustion(ci));
            for n = 1:Nfound-1
                ci = data_index(n);
                if liquid_heat_capacity_equation_num(ci) == 1
                    cp_liq_fn = liquid_heat_capacity_equation_1;
                else
                    cp_liq_fn = liquid_heat_capacity_equation_2;
                end
                if ideal_gas_heat_capacity_equation_num(ci) == 1
                    cp_ig_fn = ideal_gas_heat_capacity_equation_1;
                elseif ideal_gas_heat_capacity_equation_num(ci) == 2
                    cp_ig_fn = ideal_gas_heat_capacity_equation_2;
                else
                    cp_ig_fn = ideal_gas_heat_capacity_equation_3;
                end
                components(n) = Component( ...
                    component_name{ci}, component_formula{ci}, ...
                    Molecular_weight_data(ci), ...
                    critical_temperature_data(ci), critical_pressure_data(ci), ...
                    critical_volume_data(ci), critical_compressibility_factor(ci), ...
                    acentric_factor(ci), ...
                    vapor_pressure_equation, vapor_pressure_coefs(ci,:), ...
                    vapor_pressure_T_range(ci,:), vapor_pressure_p_range(ci,:), ...
                    dh_vaporization_equation, dh_vaporization_coefs(ci,:), ...
                    dh_vaporization_T_range(ci,:), dh_vaporization_dhv_range(ci,:), ...
                    cp_liq_fn, liquid_heat_capacity_coef(ci,:), ...
                    liquid_heat_capacity_T_range(ci,:), liquid_heat_capacity_cp_range(ci,:), ...
                    cp_ig_fn, ideal_gas_heat_capacity_coefs(ci,:), ...
                    ideal_gas_heat_capacity_T_range(ci,:), ideal_gas_heat_capacity_cp_range(ci,:), ...
                    dh_formation_ig(ci), dg_formation_ig(ci), ...
                    ds_ideal_gas(ci), dh_combustion(ci));
            end
        end
    end
end
