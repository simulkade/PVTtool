classdef Component
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
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
        function obj = Component(formula, ...
        MW, ...
        Tc, ...
        Pc, ...
        Vc, ...
        Zc, ...
        acentric_factor, ...
        Psat_eq, ...
        Psat_coefs, ...
        PsatTrange, ...
        Psatrange, ...
        dh_vap_eq, ...
        dh_vap_coefs, ...
        dh_vap_Trange, ...
        dh_vap_range, ...
        cp_liq_eq, ...
        cp_liq_coefs, ...
        cp_liq_Trange, ...
        cp_liq_range, ...
        cp_ig_eq, ...
        cp_ig_coefs, ...
        cp_ig_Trange, ...
        cp_ig_range, ...
        dhf_ig, ...
        dgf_ig, ...
        ds_ig, ...
        dh_comb)
            %Component Construct an instance of this class
            % Creates a Component object; currently it can be created from
            % a database from the Perry's handbook
            obj.name = name;
            obj.formula = formula;
            obj.MW = MW;
            obj.Tc = Tc;
            obj.Pc = Pc;
            obj.Vc = Vc;
            obj.Zc = Zc;
            obj.acentric_factor = acentric_factor;
            obj.Psat_eq = Psat_eq;
            obj.Psat_coefs = Psat_coefs;
            obj.PsatTrange = PsatTrange;
            obj.Psatrange = Psatrange;
            obj.dh_vap_eq = dh_vap_eq;
            obj.dh_vap_coefs = dh_vap_coefs;
            obj.dh_vap_Trange = dh_vap_Trange;
            obj.dh_vap_range = dh_vap_range;
            obj.cp_liq_eq = cp_liq_eq;
            obj.cp_liq_coefs = cp_liq_coefs;
            obj.cp_liq_Trange = cp_liq_Trange;
            obj.cp_liq_range = cp_liq_range;
            obj.cp_ig_eq = cp_ig_eq;
            obj.cp_ig_coefs = cp_ig_coefs;
            obj.cp_ig_Trange = cp_ig_Trange;
            obj.cp_ig_range = cp_ig_range;
            obj.dhf_ig = dhf_ig;
            obj.dgf_ig = dgf_ig;
            obj.ds_ig = ds_ig;
            obj.dh_comb = dh_comb;
        end
        
        function p_sat = vapor_pressure(obj,T)
            %vapor_pressure(T)
            %   p_sat = component.vapor_pressure(T)
            %   T is temperature in K
            %   p_sat is pressure in Pa
            p_sat = obj.Psat_eq(T, obj.Psat_coefs);
        end
        
        function dh_v = dh_vap(obj,T)
            %dh_vap(T)
            %   dh_v [J/mol] = component.dh_vap(T)
            %   T is temperature in K
            dh_v = obj.dh_vap_eq(T/obj.Tc, obj.dh_vap_coefs);
        end
        
        function cp_l = cp_liq(obj,T)
            %cp_liq(T)
            %   cp_l [J/(mol.K)] = component.cp_liq(T)
            %   liquid phase specific heat capacity
            %   T is temperature in K
            cp_l = obj.cp_liq_eq(T/obj.Tc, obj.cp_liq_coefs);
        end
        
        function cp = cp_ig(obj,T)
            %cp_ig(T)
            %   Ideal gas heat capacity
            %   cp [J/(mol.K)] = component.cp_ig(T)
            %   T is temperature in K
            cp = obj.cp_ig(T, obj.cp_ig_coefs);
        end
    end
end

