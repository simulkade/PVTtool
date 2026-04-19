function [vf, RRflag] = RachfordRiceNR(composition, ki, vapor_frac_est)
% RachfordRiceNR  Solve Rachford-Rice equation by Newton-Raphson.
%
%   [vf, RRflag] = RachfordRiceNR(composition, ki, vapor_frac_est)
%
%   Solves  f(V) = sum_i z_i*(K_i-1)/(1+V*(K_i-1)) = 0
%   for vapor fraction V, using Newton-Raphson with bounds enforcement.
%   V is constrained to the physical interval (Vmin, Vmax) determined by
%   the K-values to prevent singularities.
%
% PARAMETERS:
%   composition    - [1 x N] feed mole fractions z_i
%   ki             - [1 x N] K-values K_i = y_i/x_i
%   vapor_frac_est - initial guess for V
%
% RETURNS:
%   vf     - solved vapor fraction
%   RRflag - 1 if converged, 0 if maximum iterations reached

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

% Physical interval: V must lie in (Vmin, Vmax) to keep all denominators finite
Kmax = max(ki);
Kmin = min(ki);
if Kmax > 1
    Vmin = 1/(1-Kmax);   % most negative permissible V
else
    Vmin = -1e10;
end
if Kmin < 1
    Vmax = 1/(1-Kmin);   % largest permissible V
else
    Vmax = 1e10;
end

vapor_frac = min(max(vapor_frac_est, Vmin + 1e-10), Vmax - 1e-10);
RRflag = 0;
count  = 0;

while count < 100
    count = count + 1;
    f    = 0;
    dfdv = 0;
    for i = 1:length(composition)
        if composition(i) ~= 0
            denom = 1 + vapor_frac*(ki(i)-1);
            f    = f    + composition(i)*(ki(i)-1)/denom;
            dfdv = dfdv - composition(i)*(ki(i)-1)^2/denom^2;
        end
    end
    if abs(dfdv) < 1e-30
        break
    end
    dvf = f/dfdv;
    vapor_frac = vapor_frac - dvf;
    % Clamp back to physical interval
    vapor_frac = min(max(vapor_frac, Vmin + 1e-10), Vmax - 1e-10);
    % Convergence: relative + absolute tolerance
    if abs(dvf) < 1e-10 * (1 + abs(vapor_frac))
        RRflag = 1;
        break
    end
end
vf = vapor_frac;
