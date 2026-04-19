function [vapor_y, liquid_x, vapor_frac] = vleflash(mixture, thermo, options)
% vleflash  Isothermal Vapor-Liquid Equilibrium flash at fixed T and P.
%
%   [vapor_y, liquid_x, vapor_frac] = vleflash(mixture, thermo, options)
%
%   Solves the Rachford-Rice equation by Newton-Raphson successive
%   substitution.  Vapor fraction is constrained to [0, 1].
%   All composition vectors are row vectors [1 x N].
%
% PARAMETERS:
%   mixture    - Mixture object with temperature, pressure, mole_fraction,
%                components, and bip fields
%   thermo     - ThermoModel object specifying EOS and mixing rule
%   options    - FlashOptions object (or struct) with fields:
%                  accuracy  — convergence tolerance (default 1e-7)
%                  iteration — maximum iterations (default 100)
%
% RETURNS:
%   vapor_y    - [1 x N] vapor-phase mole fractions
%   liquid_x   - [1 x N] liquid-phase mole fractions
%   vapor_frac - scalar vapor mole fraction in [0, 1]
%
% EXAMPLE:
%   [comp, ~] = addComponents({'CH4', 'C10H22'});
%   mix = Mixture(comp, 300, 5e6);
%   mix.mole_fraction = [0.5 0.5];
%   thermo = ThermoModel();
%   opts = FlashOptions();
%   [y, x, V] = vleflash(mix, thermo, opts);
%
% SEE ALSO: vleflashnegative, stabilityTest, FlashOptions

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
    %Initialization of variables:
    vapor_frac=0.5;
    error1=1;
    error2=1;
    error3=1;
    j=0;
    while ((error1>eps1) || (error2>eps1) || (error3>eps1))
    j=j+1;
    if (j>max_itr)  %check the maximum number of itr to avoid infinite loop
      break;
    end
    [f, dfdv] = massbalfunc(composition, ki, vapor_frac);
    %New value for variable vapor_frac using Newton's method:
    vapor_frac=vapor_frac-f/dfdv;
    %physical constraint on vapor_frac
%     ss = 0.7;
%     dvf = f/dfdv;
%     while ((vapor_frac_new<0) || (vapor_frac_new>1))
%         dvf = ss*dvf;        
%         vapor_frac_new=vapor_frac-dvf;
%     end
%     vapor_frac = vapor_frac_new;
    if (vapor_frac<0)
      vapor_frac=0;
    elseif (vapor_frac>1)
      vapor_frac=1;
    end
    
    [liquid_x, vapor_y] = xy_calc(composition, vapor_frac, ki);
    
    error1=abs(sum(liquid_x)-1);
    error2=abs(sum(vapor_y)-1);
    %Another if statement based on my experience:
    if (abs(dfdv)<eps1)
        error3=f;
    else                  
        error3=abs(f/dfdv);
    end
    % normalize the mole fractions
    liquid_x = mynormalize(liquid_x);
    vapor_y = mynormalize(vapor_y);
    
    ki = kvalue(mixture, thermo, liquid_x, vapor_y);
    end
  
  % Final physical constraint for 1-phase result
  if ((vapor_frac==0) || (vapor_frac==1))
      vapor_y=composition;
      liquid_x=composition;
  end
  vapor_y = mynormalize(vapor_y);
  liquid_x = mynormalize(liquid_x);
