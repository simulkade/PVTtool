classdef Mixture
% Mixture  Multicomponent mixture at given temperature and pressure.
%
%   mix = Mixture(components, T_K, p_Pa)
%
%   Creates a Mixture from a Component array with equimolar composition
%   and a zeroed BIP object.  Adjust properties after construction:
%     mix.mole_fraction = [0.3 0.7];
%     mix.bip.EOScons   = [0 -0.05; -0.05 0];

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
        components       % 1-by-N Component array
        mole_fraction    % 1-by-N double (mole fractions, sum = 1)
        pressure         % scalar [Pa]
        temperature      % scalar [K]
        bip              % BIP object
    end

    methods
        function obj = Mixture(components, T_K, p_Pa)
            obj.components    = components;
            n                 = numel(components);
            obj.mole_fraction = ones(1, n) / n;
            obj.pressure      = p_Pa;
            obj.temperature   = T_K;
            obj.bip           = BIP(n);
        end
    end
end
