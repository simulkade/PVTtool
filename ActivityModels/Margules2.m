function [gErt, gama] = Margules2(temperature, x, component, BIP)
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


%NOT READY YET
% Margules 2 suffix activity coefficient method gives gE/RT in []
% Temperature in K
% A in J/mole  (main diagonal of the matrix should be 1)

A = [BIP.NRTLcons]+[BIP.NRTLtdep]*temperature;
alfa = BIP.NRTLalfa;
R = 8.314;  % J/mole/K
N = length(x);
mult1 = ones(N)-diag(ones(N,1));
taw = (A/R/temperature) .* mult1;
G = (exp(-alfa*taw)) .* mult1;

if (N==2)
    gErt = x(1)*x(2)*(taw(2,1)*G(2,1)/(x(1)+x(2)*G(2,1))+ ...
        taw(1,2)*G(1,2)/(x(2)+x(1)*G(1,2)));
    gama(1) = exp(x(2)^2*(taw(2,1)*(G(2,1)/(x(1)+x(2)*G(2,1)))^2+ ...
        taw(1,2)*G(1,2)/(x(2)+x(1)*G(1,2))^2));
    gama(2) = exp(x(1)^2*(taw(1,2)*(G(1,2)/(x(2)+x(1)*G(1,2)))^2+ ...
        taw(2,1)*G(2,1)/(x(1)+x(2)*G(2,1))^2));
else
    gErt = x*( x*(G.*taw)./ (x*G) )';

    part1 = x*(G.*taw)./ (x*G);
    part2 = (x.*(1./(x*G)))*(G.*taw)';
    part3 = (x.*(1./(x*G)).*(1./(x*G)).*(x*(taw.*G)))*G';
    gama = exp(part1 +part2-part3);
end