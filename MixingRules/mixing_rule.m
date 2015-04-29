function [a, b] = mixing_rule(mixture, thermo, ai, bi, s1, Q)
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


R = 8.314; %J/mol/K
bip = mixture.bip;
component = mixture.components;
mixing_rule_num = thermo.mixingrule;
activitycoef = thermo.activity_model;
temperature = mixture.temperature;
x = mixture.mole_fraction;
q1 = Q(1);
q2 = Q(2);
bipeos = [bip.EOScons]+[bip.EOStdep]*temperature;
N = length(x);
if (mixing_rule_num == 1)  %simple van der Waals mixing
    b = x*bi';
    a=0;
       for i=1:N
         for j=1:N
           a=a+x(i)*x(j)*sqrt(ai(i)*ai(j))*(1-bipeos(i,j));
         end
       end
%bipeos = [bip.EOScons]+[bip.EOStdep]*temperature;    
%a = sum(sum((x'*x).*sqrt(ai'*ai).*(1-bipeos)));   
    
    
    
elseif (mixing_rule_num == 2)  %Haron Vidal mixing rule
    [gErt, gama] = activitycoef(temperature, x, component, bip);
    b = x*bi';
    a = b*(x*(ai./bi)'-gErt*R*temperature/s1);
    
elseif (mixing_rule_num == 3)   %MHV1 mixing rule
    [gErt, gama] = activitycoef(temperature, x, component, bip);
    b = x*bi';
    %q1 = -0.593;
    %q2 = 0.0;
    alphai = ai./bi/R/temperature;  % A alpha^2 + B alpha + C = 0
    alpha = (gErt +x*log(b./bi)'+q1*x*alphai')/q1;
    a = alpha * b * R * temperature;
elseif (mixing_rule_num == 4)   %MHV2 mixing rule
    [gErt, gama] = activitycoef(temperature, x, component, bip);
    b = x*bi';
    %q1 = -0.478;
    %q2 = -0.0047;
    alphai = ai./bi/R/temperature;  % A alpha^2 + B alpha + C = 0
    C = -gErt -x*log(b./bi)'-q1*x*alphai'-q2*x*(alphai.^2)';
    B = q1;
    A = q2;
    D = B^2-4*A*C;
    alpha = (-B-sqrt(D))/(2*A);
    a = real(alpha) * b * R * temperature;
elseif (mixing_rule_num == 4)   %Wong-Sandler mixing rule
    bipeos = [bip.EOScons]+[bip.EOStdep]*temperature;
    [gErt, gama] = activitycoef(temperature, x, component, bip);
    nume1=0;
    denom1 = 0;
       for i=1:N
          denom1 = denom1 + x(i)*ai(i)/bi(i);
          for j=1:N
             nume1 = nume1 + x(i)*x(j)*0.5*(bi(i)-ai(i)/(R*temperature)+ ...
                 bi(j)-ai(j)/(R*temperature))*(1-bipeos(i,j));
             a=a+x(i)*x(j)*sqrt(ai(i)*ai(j))*(1-BIPs(i,j));
          end
       end
    b = R*temperature*nume1/(R*temperature-(denom1+gErt*R*temperature/Cstar));
    a = b*R*temperature*(gErt/Cstar+denom1/(R*temperature));
end
    
