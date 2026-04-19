function [P_bub, y_bub, conv_flag] = bubblePressure(mix, thermo, opts)
% bubblePressure  Bubble-point pressure at fixed temperature.
%
%   [P_bub, y_bub, conv_flag] = bubblePressure(mix, thermo, opts)
%
%   Finds the bubble-point pressure P_bub and incipient vapor composition
%   y_bub for the mixture feed composition at fixed T = mix.temperature.
%   Uses Wilson correlation for initialisation then successive substitution
%   of EOS K-values.  P update: P_new = P_old * sum(K_i * z_i).
%
% PARAMETERS:
%   mix        - Mixture object; mix.temperature is the fixed T [K];
%                mix.pressure is used as initial guess (overridden by Wilson)
%   thermo     - ThermoModel object specifying EOS and mixing rule
%   opts       - FlashOptions object with fields accuracy, iteration
%
% RETURNS:
%   P_bub      - bubble-point pressure [Pa]
%   y_bub      - [1 x N] incipient vapor mole fractions
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

T   = mix.temperature;
z   = mix.mole_fraction;
Tc  = [mix.components.Tc];
Pc  = [mix.components.Pc];
omega = [mix.components.acentric_factor];
eps_val = opts.accuracy;

% Wilson initial pressure estimate: P_bub = sum(z_i * Pc_i * exp(5.37*(1+wi)*(1-Tci/T)))
P = sum(z .* Pc .* exp(5.37*(1+omega).*(1-Tc/T)));
P = max(P, 1e2);  % safety floor

% Initial K and y from Wilson
K = (Pc/P) .* exp(5.37*(1+omega).*(1-Tc/T));
y = K.*z;
S = sum(y);
y = y / S;

conv_flag = 0;
for iter = 1:opts.iteration
    mix1 = mix;
    mix1.pressure = P;
    K_new = kvalue(mix1, thermo, z, y);
    S = sum(K_new .* z);
    % Update pressure: drives sum(K*z) → 1
    P = P * S;
    P = max(P, 1e2);
    % Update incipient vapor
    y = K_new.*z / S;
    y = y / sum(y);
    % Check convergence on both K stability and saturation condition
    err_K = max(abs(K_new - K) ./ max(abs(K), 1e-10));
    if err_K < eps_val && abs(S-1) < eps_val
        conv_flag = 1;
        break
    end
    K = K_new;
end
P_bub = P;
y_bub = y;
