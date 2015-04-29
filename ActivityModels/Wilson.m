function [gErt, gama] = Wilson(temperature, x, component, BIP)
%NOT READY YET
% Margules 2 suffix activity coefficient method gives gE/RT in []
% Temperature in K
% A in J/mole  (main diagonal of the matrix should be 1)
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


A = [BIP.Wilsoncons]+[BIP.Wilsontdep]*temperature;
R = 8.314;  % J/mole/K
N = length(x);

if (N==2)
    gErt = -x(1)*log(x(1)+A(1,2)*x(2))-x(2)*log(x(2)+A(2,1)*x(1));
    gama(1) = exp(-log(x(1)+A(1,2)*x(2))+x(2)*(A(1,2)/(x(1)+A(1,2)*x(2)) ...
        -A(2,1)/(A(2,1)*x(1)+x(2))));
    gama(2) = exp(-log(x(2)+A(2,1)*x(1))-x(1)*(A(1,2)/(x(1)+A(1,2)*x(2)) ...
        -A(2,1)/(A(2,1)*x(1)+x(2))));
else
    gErt = x*log(x*A')';

    part1 = -log(x*A')+1;
    part2 = -(x.*(1./(x*A')))*A;
    gama = exp(part1+part2);
end