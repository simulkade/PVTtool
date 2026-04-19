classdef ThermoModel
% ThermoModel  Thermodynamic model configuration for flash calculations.
%
%   thermo = ThermoModel()
%
%   Default: Peng-Robinson EOS, NRTL activity model, van der Waals mixing.
%   Override properties as needed:
%     thermo.EOS            = @PREOS   % or @SRKEOS, @PR78EOS
%     thermo.activity_model = @NRTL    % or @Wilson, @UNIQUAC, @Margules2
%     thermo.mixingrule     = 1        % 1=vdW, 2=HV, 3=MHV1, 4=MHV2
%     thermo.phase          = 1        % 1=liquid, 2=vapor
%     thermo.fugacity_switch = 1       % 1=compute fugacity, 0=skip

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

    properties
        EOS             = @PREOS
        activity_model  = @NRTL
        mixingrule      = 1
        phase           = 1
        fugacity_switch = 1
    end

    methods
        function obj = ThermoModel()
        end
    end
end
