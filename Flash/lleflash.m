function [liquid2_y, liquid1_x, liquid2_frac] = lleflash(mixture, thermo, options)
% lleflash  Liquid-Liquid Equilibrium flash at fixed T and P.
%
%   [liquid2_y, liquid1_x, liquid2_frac] = lleflash(mixture, thermo, options)
%
%   Uses successive substitution on the Rachford-Rice equation with LLE
%   K-values (Ki = fi_L1 / fi_L2, both phases evaluated as liquid).
%
%   IMPORTANT: all composition vectors are row vectors [1 x N].
%
% SEE ALSO: vleflash, kvalueLLE, stabilityLLETest

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

    eps1    = options.accuracy;
    max_itr = options.iteration;

    % Initial K-value estimate from Wilson correlation
    ki = kval_estimate(mixture);
    composition = mixture.mole_fraction;

    % Initialise variables
    liquid2_frac = 0.5;
    error1 = 1;
    error2 = 1;
    error3 = 1;
    j = 0;

    while ((error1 > eps1) || (error2 > eps1) || (error3 > eps1))
        j = j + 1;
        if (j > max_itr)
            break;
        end

        [f, dfdv] = massbalfunc(composition, ki, liquid2_frac);

        % Newton step on Rachford-Rice
        liquid2_frac = liquid2_frac - f/dfdv;

        % Physical bounds
        if (liquid2_frac < 0)
            liquid2_frac = 0;
        elseif (liquid2_frac > 1)
            liquid2_frac = 1;
        end

        [liquid1_x, liquid2_y] = xy_calc(composition, liquid2_frac, ki);

        error1 = abs(sum(liquid1_x) - 1);
        error2 = abs(sum(liquid2_y) - 1);
        if (abs(dfdv) < eps1)
            error3 = f;
        else
            error3 = abs(f/dfdv);
        end

        liquid1_x = mynormalize(liquid1_x);
        liquid2_y = mynormalize(liquid2_y);

        % Update K-values using LLE fugacity ratios (both phases liquid)
        ki = kvalueLLE(mixture, thermo, liquid1_x, liquid2_y);
    end

    % Single-phase result
    if ((liquid2_frac == 0) || (liquid2_frac == 1))
        liquid2_y = composition;
        liquid1_x = composition;
    end

    liquid2_y = mynormalize(liquid2_y);
    liquid1_x = mynormalize(liquid1_x);
