function [gErt, gama] = UNIQUAC(temperature, x, component, BIP)
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


N = length(x);
R = 8.314;
z = 10;  %coordination number
r= [component.uniquacR];
q= [component.uniquacQ];
A = [BIP.UNIQUACcons]+[BIP.UNIQUACtdep]*temperature+[BIP.UNIQUACtdep2]*temperature*temperature+ ...
    [BIP.UNIQUACtdepm1]/temperature+[BIP.UNIQUACtdeplog]*log(temperature);
mult1 = ones(N)-diag(ones(N,1));
fay = (x.*r)/(x*r');
teta = (x.*q)/(x*q');
taw = exp(-A/R/temperature).*mult1;
l = z/2*(r-q)-(r-1);
if (N==2)
    lngama = zeros(1,N);
    gEcrt = x(1)*log(fay(1)/x(1))+x(2)*log(fay(2)/x(2))+ ...
        z/2*(q(1)*x(1)*log(teta(1)/fay(1))+q(2)*x(2)*log(teta(2)/fay(2)));
    gErrt = -q(1)*x(1)*log(teta(1)+teta(2)*taw(2,1))- ...
             q(2)*x(2)*log(teta(2)+teta(1)*taw(1,2));
    gErt = gEcrt + gErrt;
    for i=1:N
        for j=1:N
            if (i~=j)
                lngama(i) = log(fay(i)/x(i))+z/2*q(i)*log(teta(i)/fay(i)) ...
                    + fay(j)*(l(i)-r(i)/r(j)*l(j))-q(i)*log(teta(i)+teta(j)*taw(j,i)) ...
                    + teta(j)*q(i)*(taw(j,i)/(teta(i)+teta(j)*taw(j,i))- ...
                    taw(i,j)/(teta(j)+teta(i)*taw(i,j)));
            end
        end
    end
    gama = exp(lngama);
else
    gErt = x*log(fay./x)' + z/2*x*(q.*log(teta./fay))'-x*(q.*log(teta*taw'))';
    lngama = log(fay./x) + z/2*q.*log(teta./fay) + l - fay./x*(x*l') ...
        - q.*(teta*taw') + q - q.*((1./(teta*taw))*(teta*taw')');
    gama = exp(lngama);
end
