function [stability_flag, SL, SV] = stabilityTest(mixture, thermo, varargin)
% Michelsen stability test; I have used the algorithm described here:
% https://www.e-education.psu.edu/png520/m17_p7.html
% 
% SYNOPSIS:
%   
% 
% PARAMETERS:
%   
% 
% RETURNS:
%   
% 
% EXAMPLE:
% 
% SEE ALSO:
%     

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

% extract the data from flash options
if nargin>2
    trivial_eps = varargin{1}.trivialSolutionMaxError;
    convergence_eps = varargin{1}.convergenceMaxError;
    max_itr = varargin{1}.maxIteration;
else
    % default values
    trivial_eps = 1e-5;
    convergence_eps = 1e-10;
    max_itr = 50;
end

% extract EOS function
eosf = thermo.EOS;

% switch on the fugacity calculation in thermo structure
thermo.fugacity_switch=1; % switch on

% initialize pseudo second phases
gasmixture = mixture;
liquidmixture = mixture;
gasthermo = thermo;
gasthermo.phase = 1;
liquidthermo = thermo;
liquidthermo.phase = 1;


% extract the total composition and pressure
composition = mixture.mole_fraction;
p = mixture.pressure;

% --------------------------- FIRST TEST ----------------------------------
% calculate the fugacity of the mixture, assuming it is a liquid
[~,~,fug_coef,~]=eosf(mixture, liquidthermo);
mixfug = fug_coef.*composition*p; %[Pa]

% Initial estimate for k-values using an empirical equation
ki = kval_estimate(mixture);

% assign large number to error values to begin the loop
conv_error = 1+convergence_eps;
triv_error = 1+trivial_eps;
j = 0;
while (conv_error>convergence_eps) && (triv_error>trivial_eps) && (j<max_itr)
    j = j+1; % loop counter
    % create a vapor-like second phase
    Yi = composition.*ki;
    SV = sum(Yi);

    % normalize the vapor-like mole fractions
    yi = Yi/SV;

    % calculate the fugacity of the vapor-like phase using the thermo structure
    gasmixture.mole_fraction = yi;
    [~,~,fug_coef,~]=eosf(gasmixture, gasthermo);
    gasfug = fug_coef.*yi*p; %[Pa]

    % correct K-values
    Ri = mixfug./gasfug*(1/SV);
    ki = ki.*Ri;

    % calculate the convergence and trivial solution error values
    conv_error = sum((Ri-1).^2);
    triv_error = sum(log(ki(ki>0)).^2);
end

% analyze the first test results
if (conv_error>convergence_eps)
    stability_flag(1) = 1; % converged to trivial solution
elseif (triv_error>trivial_eps)
    stability_flag(1) = 2; % converged
elseif (j>=max_itr)
    stability_flag(1) = 3; % maximum iteration reached
end

% --------------------------- SECOND TEST ---------------------------------
% calculate the fugacity of the mixture, assuming it is a gas
[~,~,fug_coef,~]=eosf(mixture, gasthermo);
mixfug = fug_coef.*composition*p; %[Pa]

% Initial estimate for k-values using an empirical equation
ki = kval_estimate(mixture);

% assign large number to error values to begin the loop
conv_error = 1+convergence_eps;
triv_error = 1+trivial_eps;
j = 0;
while (conv_error>convergence_eps) && (triv_error>trivial_eps) && (j<max_itr)
    j = j+1; % loop counter
    % create a liquid-like second phase
    Xi = composition./ki;
    SL = sum(Xi);

    % normalize the liquid-like mole fractions
    xi = Xi/SL;

    % calculate the fugacity of the liquid-like phase using the thermo structure
    liquidmixture.mole_fraction = xi;
    [~,~,fug_coef,~]=eosf(liquidmixture, liquidthermo);
    liquidfug = fug_coef.*xi*p; %[Pa]

    % correct K-values
    Ri = liquidfug./mixfug*SL;
    ki = ki.*Ri;

    % calculate the convergence and trivial solution error values
    conv_error = sum((Ri-1).^2);
    triv_error = sum(log(ki(ki>0)).^2);
end

% analyze the first test results
if (conv_error>convergence_eps)
    stability_flag(2) = 1; % converged to trivial solution
elseif (triv_error>trivial_eps)
    stability_flag(2) = 2; % converged
elseif (j>=max_itr)
    stability_flag(2) = 3; % maximum iteration reached
end
