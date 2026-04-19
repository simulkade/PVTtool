function [liquid_z, vapor_z, fugacity, HR, props] = PREOS(mixture, thermo)
% PREOS  Peng-Robinson Equation of State (1976).
%
%   [liquid_z, vapor_z, fugacity, HR] = PREOS(mixture, thermo)
%   [liquid_z, vapor_z, fugacity, HR, props] = PREOS(mixture, thermo)
%
%   Computes compressibility factors, fugacity coefficients, residual
%   enthalpy, and (optionally) a full set of residual properties for a
%   mixture at the conditions stored in MIXTURE.
%   Supports van der Waals (rule 1), Huron-Vidal (rule 2), MHV1 (rule 3),
%   and MHV2 (rule 4) mixing rules.
%
% PARAMETERS:
%   mixture  - Mixture object; relevant fields:
%              .temperature  [K]
%              .pressure     [Pa]
%              .mole_fraction [1 x N]
%              .components   [1 x N Component array]  (Tc, Pc, acentric_factor)
%              .bip          BIP object (EOScons, EOStdep for rule 1)
%   thermo   - ThermoModel object; relevant fields:
%              .mixingrule       1=vdW, 2=HV, 3=MHV1, 4=MHV2
%              .activity_model   function handle (used for rules 2-4)
%              .phase            1=liquid (use Zl), 2=vapor (use Zv)
%              .fugacity_switch  1=compute fugacity, 0=skip
%
% RETURNS:
%   liquid_z  - liquid compressibility factor (smallest real root > B)
%   vapor_z   - vapor compressibility factor (largest real root > B)
%   fugacity  - [1 x N] fugacity coefficients  (zero if fugacity_switch==0)
%   HR        - residual molar enthalpy [J/mol]
%   props     - struct with fields HR, SR, GR, VR, Cp_R, Cv_R  (optional 5th output)
%               SR [J/(mol·K)], GR [J/mol], VR [m³/mol], Cp_R/Cv_R [J/(mol·K)]
%
% SEE ALSO: SRKEOS, PR78EOS, select_z_roots, ThermoModel

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

mixing_rule_num = thermo.mixingrule;
activityfun = thermo.activity_model;
phase1 = thermo.phase;
fug_need = thermo.fugacity_switch;
critical_pres = [mixture.components.Pc]; %[Pa]
critical_temp = [mixture.components.Tc]; %[K]
acentric_fact = [mixture.components.acentric_factor]; %[-]
BIP = mixture.bip;
x = mixture.mole_fraction;
p = mixture.pressure; %[Pa]
T = mixture.temperature;
N = length(critical_temp);
fugacity = zeros(1,N);
R = 8.314;
s1 = 0.623225;  %Huron Vidal constant for PR

bi  = 0.077796*R*critical_temp./critical_pres;
aci = 0.457235*(R*critical_temp).^2 ./critical_pres;
mi  = 0.37646+(1.54226-0.26992*acentric_fact).*acentric_fact;
Tr  = T./critical_temp;
alfai = 1+mi.*(1-sqrt(Tr));
alfa  = alfai.^2;
ai    = aci .* alfa;

Q = (mixing_rule_num==3)*[-0.53 0]+(mixing_rule_num==4)*[-0.4347 -0.003654];
[a, b] = mixing_rule(mixture, thermo, ai, bi, s1, Q);

A_coef = a*p/(R*T)^2;
B_coef = b*p/(R*T);

poly_coef = [1 -1+B_coef A_coef-B_coef*(2.0+3*B_coef) ...
    -B_coef*(A_coef-B_coef*(1+B_coef))];

z_root = roots(poly_coef);
[liquid_z, vapor_z] = select_z_roots(z_root, B_coef);

if (phase1 == 1)
    zz = liquid_z;
else
    zz = vapor_z;
end

if (fug_need == 1)
    if (mixing_rule_num == 1)
        part1 = bi/b*(zz-1) - log(zz-b*p/(R*T));
        part2 = x*(sqrt(ai'*ai).*(1-[BIP.EOScons]-[BIP.EOStdep]*T))';
        part3 = A_coef/(2.828*B_coef)*(bi/b-2/a*part2) ...
            *log((zz+2.414*b*p/(R*T))/(zz-0.414*b*p/(R*T)));
        fugacity = exp(part1+part3);
    elseif (mixing_rule_num == 2)
        [~, gama] = activityfun(T, x, mixture.components, BIP);
        part1 = bi/b*(zz-1) - log(zz-b*p/(R*T));
        part3 = -1/(2*sqrt(2))*(ai./bi/R/T - ...
                   log(gama)/0.623225)* ...
                   log((zz+(1+sqrt(2))*B_coef)/(zz+(1-sqrt(2))*B_coef));
        fugacity = exp(part1+part3);
    elseif (mixing_rule_num == 3)
        [~, gama] = activityfun(T, x, mixture.components, BIP);
        q1 = -0.53;
        logfi = bi/b*(zz-1) - log(zz-B_coef) - 1/(2*sqrt(2))*(ai./(bi*R*T) ...
            + log(gama)/q1+log(b./bi)/q1+(bi/b-1)/q1)*log((zz+(1+sqrt(2))*B_coef)/ ...
            (zz+(1-sqrt(2))*B_coef));
        fugacity = exp(logfi);
    elseif (mixing_rule_num == 4)
        [~, gama] = activityfun(T, x, mixture.components, BIP);
        q1 = -0.4347;
        q2 = -0.003654;
        alphai = ai./(bi*R*T);
        alpha  = a/(b*R*T);
        logfi  = bi/b*(zz-1) - log(zz-B_coef) - ...
            1/(2*sqrt(2))*(q1*alphai+q2*(alpha^2+alphai.^2)+log(gama) ...
            +log(b./bi)+bi/b-1)/(q1+2*q2*alpha)* ...
            log((zz+(1+sqrt(2))*B_coef)/(zz+(1-sqrt(2))*B_coef));
        fugacity = exp(logfi);
    end
end

%% Residual properties
% da/dT (kij=0 approximation, consistent with standard PR departure functions)
% d(sqrt(a_i))/dT = -sqrt(aci_i)*mi/(2*sqrt(Tci*T))
% da/dT (vdW, kij=0) = 2*(sum xi*sqrt(ai)) * sum(xi * d(sqrt(ai))/dT)
dadT = -(sum(x.*sqrt(ai))) .* sum(x.*sqrt(aci).*mi./sqrt(critical_temp)) ./ sqrt(T);

ln_term = log((zz+B_coef*(1+sqrt(2))) / (zz+B_coef*(1-sqrt(2))));
HR = R*T*(zz-1) + (T*dadT - a) / (b*2*sqrt(2)) * ln_term;

if nargout >= 5
    SR  = R*log(zz - B_coef) + dadT / (b*2*sqrt(2)) * ln_term;
    GR  = HR - T*SR;
    VR  = R*T*(zz-1) / p;

    % Second T-derivative of a for Cv_R and Cp_R
    % d2(sqrt(ai))/dT2 = sqrt(aci)*mi / (4*sqrt(Tci)*T^(3/2))
    dsqrtaidT   = -sqrt(aci) .* mi ./ (2*sqrt(critical_temp*T));
    d2sqrtaidT2 =  sqrt(aci) .* mi ./ (4*sqrt(critical_temp) .* T^(3/2));
    d2adT2 = 2*(sum(x.*dsqrtaidT))^2 + 2*sum(x.*sqrt(ai))*sum(x.*d2sqrtaidT2);

    Cv_R = -T * d2adT2 / (b*2*sqrt(2)) * ln_term;

    % Cp_R via thermodynamic identity: Cp - Cv = -T*(dP/dV)^{-1}*(dP/dT)^2
    % Cp_R = Cv_R + (Cp-Cv) - R  since Cp^ig - Cv^ig = R
    V_mol  = zz*R*T / p;
    den_pr = V_mol^2 + 2*b*V_mol - b^2;          % (V+b(1+√2))(V+b(1-√2)) = V²+2bV-b²
    dPdT_V = R/(V_mol-b) - dadT/den_pr;
    dPdV_T = -R*T/(V_mol-b)^2 + 2*a*(V_mol+b)/den_pr^2;
    Cp_R   = Cv_R - T*dPdT_V^2/dPdV_T - R;

    props = struct('HR',HR, 'SR',SR, 'GR',GR, 'VR',VR, 'Cp_R',Cp_R, 'Cv_R',Cv_R);
end
