function [P_dew, x_dew, conv_flag] = dewPressure(mix, thermo, opts)
% dewPressure  Dew-point pressure at fixed temperature.
%
%   [P_dew, x_dew, conv_flag] = dewPressure(mix, thermo, opts)
%
%   Finds the dew-point pressure P_dew and incipient liquid composition
%   x_dew for the mixture feed (vapor = feed z) at fixed T = mix.temperature.
%   Uses Wilson correlation for initialisation then successive substitution
%   of EOS K-values.  P update: P_new = P_old / sum(z_i / K_i).
%
% PARAMETERS:
%   mix        - Mixture object; mix.temperature is the fixed T [K];
%                mix.pressure is used as starting neighbourhood
%   thermo     - ThermoModel object specifying EOS and mixing rule
%   opts       - FlashOptions object with fields accuracy, iteration
%
% RETURNS:
%   P_dew      - dew-point pressure [Pa]
%   x_dew      - [1 x N] incipient liquid mole fractions
%   conv_flag  - 1 if converged, 0 if maximum iterations reached

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

T     = mix.temperature;
z     = mix.mole_fraction;
Tc    = [mix.components.Tc];
Pc    = [mix.components.Pc];
omega = [mix.components.acentric_factor];
eps_val = opts.accuracy;

% Wilson initial dew pressure: P_dew = 1/sum(z_i/K_wilson_i) at P=1
% Since K_wilson = Pc/P * exp(...), sum(z/K) = P/sum(z*exp(...)*Pc)
% => P_dew_Wilson = 1 / sum(z_i / (Pc_i*exp(5.37*(1+wi)*(1-Tci/T))))
K_coef = Pc .* exp(5.37*(1+omega).*(1-Tc/T));   % K_wilson * P
P = 1 / sum(z ./ K_coef);
P = max(P, 1e2);

% Initial x and K from Wilson
K = K_coef / P;
x = (z ./ K);
x = x / sum(x);

conv_flag = 0;
for iter = 1:opts.iteration
    mix1 = mix;
    mix1.pressure = P;
    K_new = kvalue(mix1, thermo, x, z);  % liquid=x, vapor=z
    H = sum(z ./ K_new);
    % Update pressure: P_new = P_old / H drives H → 1
    P = P / H;
    P = max(P, 1e2);
    % Update incipient liquid
    x = z ./ K_new;
    x = x / sum(x);

    err_K = max(abs(K_new - K) ./ max(abs(K), 1e-10));
    if err_K < eps_val && abs(H-1) < eps_val
        conv_flag = 1;
        break
    end
    K = K_new;
end
P_dew = P;
x_dew = x;
