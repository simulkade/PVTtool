function ki = kvalueLLE(mixture, thermo, liquid_x, vapor_y)
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


thermo.phase=1;
thermo.fugacity_switch=1; % switch on
eosf = thermo.EOS;
mixture1 = mixture;
mixture1.mole_fraction = liquid_x;
[zl,zv,liq_fug,hh]=eosf(mixture1, thermo);
thermo.phase=1;
mixture1.mole_fraction = vapor_y;
[zl,zv,vap_fug,hh]=eosf(mixture1, thermo);
% [zl,zv,vap_fug,hh]=EOS(critical_temp,critical_pres,vapor_y, ...
%     acentric_fact,[0 0;0 0],pressure,temperature,eos_type,1,2);
N = length(liquid_x);
    ki = zeros(1, N);
    for i=1:N
      if (vap_fug(i)~=0)
        ki(i)=liq_fug(i)/vap_fug(i);
      end
    end