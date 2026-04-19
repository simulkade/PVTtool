function [liquid_z, vapor_z, fugacity, HR, props] = SRKEOS(mixture, thermo)
% SRKEOS  Soave-Redlich-Kwong Equation of State (1972).
%
%   [liquid_z, vapor_z, fugacity, HR] = SRKEOS(mixture, thermo)
%   [liquid_z, vapor_z, fugacity, HR, props] = SRKEOS(mixture, thermo)
%
%   Computes compressibility factors, fugacity coefficients, and residual
%   properties using the SRK EOS.  Supports van der Waals (rule 1),
%   Huron-Vidal (rule 2), MHV1 (rule 3), and MHV2 (rule 4) mixing rules.
%
% PARAMETERS:
%   mixture  - Mixture object; relevant fields:
%              .temperature, .pressure, .mole_fraction, .components, .bip
%   thermo   - ThermoModel object; relevant fields:
%              .mixingrule, .activity_model, .phase, .fugacity_switch
%
% RETURNS:
%   liquid_z  - liquid compressibility factor (smallest real root > B)
%   vapor_z   - vapor compressibility factor (largest real root > B)
%   fugacity  - [1 x N] fugacity coefficients  (zero if fugacity_switch==0)
%   HR        - residual molar enthalpy [J/mol]
%   props     - struct with fields HR, SR, GR, VR, Cp_R, Cv_R  (optional 5th output)
%
% SEE ALSO: PREOS, PR78EOS, select_z_roots, ThermoModel

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
s1 = log(2);  % Huron-Vidal constant for SRK

bi  = 0.08664*R*critical_temp./critical_pres;
aci = 0.42748*(R*critical_temp).^2 ./ critical_pres;
mi  = 0.48+(1.574-0.176*acentric_fact).*acentric_fact;
Tr  = T./critical_temp;
alfai = 1+mi.*(1-sqrt(Tr));
alfa  = alfai.^2;
ai    = aci.*alfa;

Q = (mixing_rule_num==3)*[-0.593 0]+(mixing_rule_num==4)*[-0.478 -0.0047];
[a, b] = mixing_rule(mixture, thermo, ai, bi, s1, Q);

A_coef = a*p/(R*T)^2;
B_coef = b*p/(R*T);

poly_coef = [1 -1 A_coef-B_coef*(1+B_coef) -A_coef*B_coef];

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
        part3 = A_coef/B_coef*(bi/b-2/a*part2)*log(1+b*p/(R*T*zz));
        fugacity = exp(part1+part3);
    elseif (mixing_rule_num == 2)
        alphai = ai./(bi*R*T);
        alpha  = a/(b*R*T);
        [~, gama] = activityfun(T, x, mixture.components, BIP);
        alphabari = alphai - log(gama)/s1;
        logfi = log(1/(zz-B_coef)) + ((p/(R*T))/(zz-B_coef) ...
            -(p*alpha/(R*T))/(zz+B_coef))*bi - ...
            alphabari*log((zz+B_coef)/zz);
        fugacity = exp(logfi);
    elseif (mixing_rule_num == 3)
        q1 = -0.593;
        q2 = 0.0;
        alphai = ai./(bi*R*T);
        alpha  = a/(b*R*T);
        [~, gama] = activityfun(T, x, mixture.components, BIP);
        alphabari = 1/(q1+2*alpha*q2)*(q1*alphai+q2*(alpha^2+alphai.^2)+ ...
            log(gama)+log(b./bi)+bi/b-1);
        logfi = log(1/(zz-B_coef)) + ((p/(R*T))/(zz-B_coef) ...
            -(p*alpha/(R*T))/(zz+B_coef))*bi - ...
            alphabari*log((zz+B_coef)/zz);
        fugacity = exp(logfi);
    elseif (mixing_rule_num == 4)
        q1 = -0.478;
        q2 = -0.0047;
        alphai = ai./(bi*R*T);
        alpha  = a/(b*R*T);
        [~, gama] = activityfun(T, x, mixture.components, BIP);
        alphabari = 1/(q1+2*alpha*q2)*(q1*alphai+q2*(alpha^2+alphai.^2)+ ...
            log(gama)+log(b./bi)+bi/b-1);
        logfi = log(1/(zz-B_coef)) + ((p/(R*T))/(zz-B_coef) ...
            -(p*alpha/(R*T))/(zz+B_coef))*bi - ...
            alphabari*log((zz+B_coef)/zz);
        fugacity = exp(logfi);
    end
end

%% Residual properties (SRK)
% da/dT (kij=0 approximation)
% SRK integral: ln((Z+B)/Z) = ln(1 + B/Z)
dadT = -(sum(x.*sqrt(ai))) .* sum(x.*sqrt(aci).*mi./sqrt(critical_temp)) ./ sqrt(T);

ln_term = log((zz + B_coef) / zz);   % ln(1 + B/Z) for SRK
HR = R*T*(zz-1) + (T*dadT - a) / b * ln_term;

if nargout >= 5
    SR  = R*log(zz - B_coef) + dadT / b * ln_term;
    GR  = HR - T*SR;
    VR  = R*T*(zz-1) / p;

    dsqrtaidT   = -sqrt(aci) .* mi ./ (2*sqrt(critical_temp*T));
    d2sqrtaidT2 =  sqrt(aci) .* mi ./ (4*sqrt(critical_temp) .* T^(3/2));
    d2adT2 = 2*(sum(x.*dsqrtaidT))^2 + 2*sum(x.*sqrt(ai))*sum(x.*d2sqrtaidT2);

    Cv_R = -T * d2adT2 / b * ln_term;

    V_mol  = zz*R*T / p;
    den_srk = V_mol * (V_mol + b);               % V*(V+b) denominator for SRK
    dPdT_V  = R/(V_mol-b) - dadT/den_srk;
    dPdV_T  = -R*T/(V_mol-b)^2 + a*(2*V_mol+b)/den_srk^2;
    Cp_R    = Cv_R - T*dPdT_V^2/dPdV_T - R;

    props = struct('HR',HR, 'SR',SR, 'GR',GR, 'VR',VR, 'Cp_R',Cp_R, 'Cv_R',Cv_R);
end
