function [vapor_y, liquid_x, vapor_frac]=vleflashnegative(mixture, thermo, options)
%IMPORTANT: every variable should be in the form of row vectors [1 x N]
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

    eps1 = options.accuracy; 
    max_itr = options.iteration;
    %Initial estimate for k-values using an empirical equation
    ki = kval_estimate(mixture);
    composition = mixture.mole_fraction;
    [ki_max, max_index] = max(ki);
    % internal step: Correct the maximum ki
    if (ki_max<1)
            ki(max_index) = 1 + ki_max;
    end
    [ki_min, min_index] = min(ki);
    if (ki_min>1)
            ki(min_index) = 0.1;
    end
    vapor_frac = 0.5;%*(vap_frac_min+vap_frac_max);
    while(1)
        % Step 2a: find min and max K-value
        ki_min = min(ki);
        ki_max = max(ki);
        % Step 2b: find min and max vapor fractions
        vap_frac_min = 1/(1-ki_max);
        vap_frac_max = 1/(1-ki_min);
        % Step 3: solve Rachford-Rice equation
        % Newton-Raphson
%         NRflag =1;
%         [f, dfdv] = massbalfunc(composition, ki, vapor_frac);
%         vapor_frac=vapor_frac-f/dfdv;
         [vapor_frac, NRflag] = RachfordRiceNR(composition, ki, 0.5*(vap_frac_min+vap_frac_max));
        % Bi-section
        if ((vapor_frac<vap_frac_min) || (vapor_frac>vap_frac_max) || (NRflag==0))
            f=@(vapor_frac)(sum(composition.*(ki-1)./(1+vapor_frac*(ki-1))));
            vapor_frac = bisection(f, ...
                vap_frac_min+eps, ...
                vap_frac_max-eps);
%             vapor_frac = fzero(f, ...
%                 [vap_frac_min+0.02*(vap_frac_max-vap_frac_min) vap_frac_max-0.02*(vap_frac_max-vap_frac_min)]);
%               vapor_frac = 0.5*(vap_frac_min+vap_frac_max);
        end  
        % Step 4: calculate x and y and normilize
        [liquid_x, vapor_y] = xy_calc(composition, vapor_frac, ki);
        liquid_x = mynormalize(liquid_x);
        vapor_y = mynormalize(vapor_y);
        % Step 5: calculate fugacity values for liquid and gas phases
        [liq_fug, vap_fug] = fugacity(mixture, thermo, liquid_x, vapor_y);

        % Step 6: check the equality of fugacities
        fug_index = (vap_fug~=0);
        error1 = sum(abs(liq_fug(fug_index)./vap_fug(fug_index)-1));
        if (error1<eps1)
            break
        end
        if (sum(log(ki).^2)<1e-4)
            if (vapor_frac<0)
               vapor_frac=0;
            elseif (vapor_frac>1)
               vapor_frac=1;
            end
            liquid_x = composition;
            vapor_y = composition;
            break
        end

        % Step 7: update K-values
%         ki = ki.*liq_fug./vap_fug;
        ki = kvalue(mixture, thermo, liquid_x, vapor_y);
        
    end
    
    
%     vapor_frac = vapor_frac+eps;
    %Initialization of variables:
%     vapor_frac=0.5;
%     error1=1;
%     error2=1;
%     error3=1;
%     j=0;
%     while ((error1>eps1) || (error2>eps1) || (error3>eps1))
%     j=j+1;
%     if (j>max_itr)  %check the maximum number of itr to avoid infinite loop
%       break;
%     end
%     [f, dfdv] = massbalfunc(composition, ki, vapor_frac);
%     %New value for variable vapor_frac using Newton's method:
%     vapor_frac=vapor_frac-f/dfdv;
%     %physical constraint on vapor_frac
% %     ss = 0.7;
% %     dvf = f/dfdv;
% %     while ((vapor_frac_new<0) || (vapor_frac_new>1))
% %         dvf = ss*dvf;        
% %         vapor_frac_new=vapor_frac-dvf;
% %     end
% %     vapor_frac = vapor_frac_new;
%     if (vapor_frac<0)
%       vapor_frac=0;
%     elseif (vapor_frac>1)
%       vapor_frac=1;
%     end
%     
%     [liquid_x, vapor_y] = xy_calc(composition, vapor_frac, ki);
%     
%     error1=abs(sum(liquid_x)-1);
%     error2=abs(sum(vapor_y)-1);
%     %Another if statement based on my experience:
%     if (abs(dfdv)<eps1)
%         error3=f;
%     else                  
%         error3=abs(f/dfdv);
%     end
%     % normalize the mole fractions
%     liquid_x = mynormalize(liquid_x);
%     vapor_y = mynormalize(vapor_y);
%     
%     ki = kvalue(component, BIP, liquid_x, vapor_y, pressure, ...
%         temperature, eosf, mixrule, activityfun);
%     end
%   
%   % Final physical constraint for 1-phase result
%   if ((vapor_frac==0) || (vapor_frac==1))
%       vapor_y=composition;
%       liquid_x=composition;
%   end
%   vapor_y = mynormalize(vapor_y);
%   liquid_x = mynormalize(liquid_x);
