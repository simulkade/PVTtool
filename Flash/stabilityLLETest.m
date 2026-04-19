function [stability_flag, SL, SV, result] = stabilityLLETest(mixture, thermo, varargin)
% stabilityLLETest  Michelsen two-sided stability test for LLE.
%
%   [stability_flag, SL, SV] = stabilityLLETest(mixture, thermo)
%   [stability_flag, SL, SV] = stabilityLLETest(mixture, thermo, options)
%   [stability_flag, SL, SV, result] = stabilityLLETest(...)
%
%   Tests whether a single liquid phase at the conditions stored in
%   MIXTURE will spontaneously split into two liquid phases using
%   Michelsen's successive-substitution algorithm.  Both trial phases are
%   evaluated as liquid (phase = 1).
%
%   When called with no output arguments the human-readable result is
%   printed to the console.
%
% PARAMETERS:
%   mixture  - Mixture object (temperature, pressure, mole fractions, bip)
%   thermo   - ThermoModel object
%   options  - (optional) FlashOptions object; relevant fields:
%              trivialSolutionMaxError  (default 1e-5)
%              convergenceMaxError      (default 1e-10)
%              maxIteration             (default 50)
%
% RETURNS:
%   stability_flag - 1-by-2 integer array; one value per test direction:
%                    1 = trivial solution (stable in this direction)
%                    2 = non-trivial convergence (unstable, split expected)
%                    3 = maximum iterations reached (inconclusive)
%   SL             - saturation of second liquid phase from Test 2
%   SV             - saturation of second liquid phase from Test 1
%   result         - struct with fields:
%                    .overall        'stable'|'unstable'|'inconclusive'
%                    .message        human-readable summary string
%                    .test1_message  string for first test direction
%                    .test2_message  string for second test direction
%                    .SL             same as SL
%                    .SV             same as SV
%
% EXAMPLE:
%   [component, ~] = addComponents({'H2O', 'C10H22'});
%   mix = Mixture(component, 300, 100e5);
%   mix.mole_fraction = [0.5 0.5];
%   thermo = ThermoModel();
%   [flag, SL, SV, res] = stabilityLLETest(mix, thermo);
%   disp(res.message)
%
% SEE ALSO: stabilityTest, lleflash, FlashOptions

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

if nargin > 2
    trivial_eps    = varargin{1}.trivialSolutionMaxError;
    convergence_eps = varargin{1}.convergenceMaxError;
    max_itr        = varargin{1}.maxIteration;
else
    trivial_eps    = 1e-5;
    convergence_eps = 1e-10;
    max_itr        = 50;
end

eosf = thermo.EOS;
thermo.fugacity_switch = 1;

phase1mixture = mixture;
phase2mixture = mixture;
liquidthermo  = thermo;  liquidthermo.phase = 1;

composition = mixture.mole_fraction;
p           = mixture.pressure;

stability_flag = zeros(1, 2);

% -------------------------------------------------------------------------
% TEST 1: can a second liquid form? (feed as liquid, trial as liquid)
% -------------------------------------------------------------------------
[~, ~, fug_coef, ~] = eosf(mixture, liquidthermo);
mixfug = fug_coef .* composition * p;

ki = kval_estimate(mixture);
conv_error = 1 + convergence_eps;
triv_error = 1 + trivial_eps;
j = 0;
while (conv_error > convergence_eps) && (triv_error > trivial_eps) && (j < max_itr)
    j = j + 1;
    Yi = composition .* ki;
    SV = sum(Yi);
    yi = Yi / SV;
    phase1mixture.mole_fraction = yi;
    [~, ~, fug_coef, ~] = eosf(phase1mixture, liquidthermo);
    phase1fug = fug_coef .* yi * p;
    Ri = mixfug ./ phase1fug * (1/SV);
    ki = ki .* Ri;
    conv_error = sum((Ri - 1).^2);
    triv_error = sum(log(ki(ki > 0)).^2);
end

if triv_error <= trivial_eps
    stability_flag(1) = 1;
    t1msg = 'Test 1 (phase A): trivial solution — no second liquid phase forms.';
elseif conv_error <= convergence_eps
    stability_flag(1) = 2;
    t1msg = 'Test 1 (phase A): non-trivial convergence — second liquid phase can form.';
else
    stability_flag(1) = 3;
    t1msg = 'Test 1 (phase A): maximum iterations reached — inconclusive.';
end

% -------------------------------------------------------------------------
% TEST 2: can a second liquid form from the other direction?
% -------------------------------------------------------------------------
[~, ~, fug_coef, ~] = eosf(mixture, liquidthermo);
mixfug = fug_coef .* composition * p;

ki = kval_estimate(mixture);
conv_error = 1 + convergence_eps;
triv_error = 1 + trivial_eps;
j = 0;
while (conv_error > convergence_eps) && (triv_error > trivial_eps) && (j < max_itr)
    j = j + 1;
    Xi = composition ./ ki;
    SL = sum(Xi);
    xi = Xi / SL;
    phase2mixture.mole_fraction = xi;
    [~, ~, fug_coef, ~] = eosf(phase2mixture, liquidthermo);
    phase2fug = fug_coef .* xi * p;
    Ri = phase2fug ./ mixfug * SL;
    ki = ki .* Ri;
    conv_error = sum((Ri - 1).^2);
    triv_error = sum(log(ki(ki > 0)).^2);
end

if triv_error <= trivial_eps
    stability_flag(2) = 1;
    t2msg = 'Test 2 (phase B): trivial solution — no second liquid phase forms.';
elseif conv_error <= convergence_eps
    stability_flag(2) = 2;
    t2msg = 'Test 2 (phase B): non-trivial convergence — second liquid phase can form.';
else
    stability_flag(2) = 3;
    t2msg = 'Test 2 (phase B): maximum iterations reached — inconclusive.';
end

% -------------------------------------------------------------------------
% Build result struct
% -------------------------------------------------------------------------
result.test1_message = t1msg;
result.test2_message = t2msg;
result.SL = SL;
result.SV = SV;

if any(stability_flag == 2)
    result.overall = 'unstable';
    result.message = sprintf( ...
        'Mixture is UNSTABLE (LLE). Liquid-liquid split expected at T=%.2f K, P=%.4g Pa.\n  %s\n  %s', ...
        mixture.temperature, mixture.pressure, t1msg, t2msg);
elseif all(stability_flag == 1)
    result.overall = 'stable';
    result.message = sprintf( ...
        'Mixture is stable (LLE). Single liquid phase at T=%.2f K, P=%.4g Pa.\n  %s\n  %s', ...
        mixture.temperature, mixture.pressure, t1msg, t2msg);
else
    result.overall = 'inconclusive';
    result.message = sprintf( ...
        'LLE stability result is inconclusive at T=%.2f K, P=%.4g Pa.\n  %s\n  %s', ...
        mixture.temperature, mixture.pressure, t1msg, t2msg);
end

if nargout == 0
    disp(result.message);
end
