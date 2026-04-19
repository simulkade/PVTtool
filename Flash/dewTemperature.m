function [T_dew, x_dew, conv_flag] = dewTemperature(mix, thermo, opts)
% dewTemperature  Dew-point temperature at fixed pressure.
%
%   [T_dew, x_dew, conv_flag] = dewTemperature(mix, thermo, opts)
%
%   Finds the dew-point temperature T_dew and incipient liquid composition
%   x_dew for the mixture feed (vapor = feed z) at fixed P = mix.pressure.
%   Uses Wilson Newton for T initialisation then successive substitution
%   of EOS K-values with a Newton T-update at each step.
%
% PARAMETERS:
%   mix        - Mixture object; mix.pressure is the fixed P [Pa];
%                mix.temperature is used as starting guess [K]
%   thermo     - ThermoModel object specifying EOS and mixing rule
%   opts       - FlashOptions object with fields accuracy, iteration
%
% RETURNS:
%   T_dew      - dew-point temperature [K]
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

P     = mix.pressure;
z     = mix.mole_fraction;
Tc    = [mix.components.Tc];
Pc    = [mix.components.Pc];
omega = [mix.components.acentric_factor];
eps_val = opts.accuracy;

% Wilson Newton iteration for initial T estimate
% At dew: H(T) = sum(z_i/K_wilson_i) = 1
% K_wilson = Pc/P * exp(5.37*(1+w)*(1-Tc/T))
% dH/dT = -sum(z/K * 5.37*(1+w)*Tc/T^2)   (always < 0)
T = mix.temperature;
for i = 1:30
    K_w  = (Pc/P) .* exp(5.37*(1+omega).*(1-Tc/T));
    H    = sum(z ./ K_w);
    dHdT = -sum(z ./ K_w .* 5.37.*(1+omega).*Tc ./ T^2);
    dT   = (H-1) / min(dHdT, -1e-20);
    T    = T - dT;
    T    = max(50, min(T, 3000));
    if abs(dT) < 0.01; break; end
end

% EOS successive substitution with T update
K = (Pc/P) .* exp(5.37*(1+omega).*(1-Tc/T));
x = z ./ K;
x = x / sum(x);

conv_flag = 0;
for iter = 1:opts.iteration
    mix1 = mix;
    mix1.temperature = T;
    K_new = kvalue(mix1, thermo, x, z);  % liquid=x, vapor=z
    H = sum(z ./ K_new);
    x = z ./ K_new;
    x = x / sum(x);

    % Newton T-update: dH/dT < 0, so step is -(H-1)/dHdT
    dHdT = -sum(z ./ K_new .* 5.37.*(1+omega).*Tc ./ T^2);
    T    = T - 0.3*(H-1) / min(dHdT, -1e-20);
    T    = max(50, min(T, 3000));

    err_K = max(abs(K_new - K) ./ max(abs(K), 1e-10));
    if err_K < eps_val && abs(H-1) < eps_val
        conv_flag = 1;
        break
    end
    K = K_new;
end
T_dew = T;
x_dew = x;
