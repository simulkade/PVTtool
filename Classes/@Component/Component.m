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
        function obj = untitled(inputArg1,inputArg2)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

