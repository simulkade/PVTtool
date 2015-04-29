function [mu_g, msg] = puregasviscosity(component, T, p, rho, method)
% T	[K]
% p	[Pa]
% rho	[kmol/m^3]
% method	'stiel-thodos'		pr<0.6, pr>0.6, hydrocarbons
%		'Reichenberg'		pr<0.6, hydrocarbons
% 		'stiel-thodos-nhc'	pr>0.6, non-hydrocarbons, polar, non-polar

% For prediction of the vapor viscosity of pure hydrocarbons at low
% pressure (below Pr of 0.6), the method of Stiel and Thodos is the
% most accurate
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


% pure component data
tc = component.Tc; % critical temperature, K
vc = component.Vc; % critical volume, m^3/kmol
rho_r = rho*vc; % reduced molar density, dimensionless
Tr = T/tc; % reduced temperature, dimensionless
M = component.MW; % kg/kmol, molecular weight
pc = component.Pc; % critical pressure, Pa
Pr = P/pc; % reduced pressure, dimensionless
switch lower(method)
	case 'stiel-thodos'
		if Tr<=1.5
			N = 0.0003400*Tr^0.94;
		else   % Tr>1.5
			N = 0.0001778*(4.58*Tr-1.67)^0.625;
	 	end
		mu_g = 4.60e-4*(N*sqrt(M)*pc^(2/3))/tc^(1/6)*1e-3; % Pa.s (1e-3 is for the conversion from mPa.s to Pa.s)
		if Pr>0.6
			mu_g = mu_g + 5.0e-8*sqrt(M)*pc^(2/3)/tc^(1/6)*(exp(1.439*rho_r)-exp(-1.11*rho_r^1.858))*1e-3; % 1e-3 like above
		end
	case 'Reichenberg'	
	case 'stiel-thodos-nhc'
end
	
