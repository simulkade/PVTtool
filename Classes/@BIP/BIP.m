classdef BIP
% BIP  Binary Interaction Parameters container.
%
%   bip = BIP(n) creates an n-component BIP object with all parameter
%   matrices initialised to zero. Set individual matrices after
%   construction to provide non-zero interaction parameters.
%
%   Properties (all n-by-n matrices unless noted):
%     EOScons, EOStdep                  - EOS kij parameters
%     NRTLcons, NRTLtdep, NRTLtdep2    - NRTL: A + B*T + C*T^2
%     NRTLtdepm1, NRTLtdeplog          - NRTL: D/T + E*ln(T) terms
%     NRTLalfa                          - NRTL non-randomness parameter
%     Wilsoncons, Wilsontdep            - Wilson parameters
%     UNIQUACcons, UNIQUACtdep         - UNIQUAC interaction parameters
%     UNIQUACR  (1-by-n)               - UNIQUAC volume parameters
%     UNIQUACQ  (1-by-n)               - UNIQUAC surface parameters

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
        EOScons
        EOStdep
        NRTLcons
        NRTLtdep
        NRTLtdep2
        NRTLtdepm1
        NRTLtdeplog
        NRTLalfa
        Wilsoncons
        Wilsontdep
        UNIQUACcons
        UNIQUACtdep
        UNIQUACR
        UNIQUACQ
    end

    methods
        function obj = BIP(n)
            obj.EOScons     = zeros(n);
            obj.EOStdep     = zeros(n);
            obj.NRTLcons    = zeros(n);
            obj.NRTLtdep    = zeros(n);
            obj.NRTLtdep2   = zeros(n);
            obj.NRTLtdepm1  = zeros(n);
            obj.NRTLtdeplog = zeros(n);
            obj.NRTLalfa    = zeros(n);
            obj.Wilsoncons  = zeros(n);
            obj.Wilsontdep  = zeros(n);
            obj.UNIQUACcons = zeros(n);
            obj.UNIQUACtdep = zeros(n);
            obj.UNIQUACR    = zeros(1, n);
            obj.UNIQUACQ    = zeros(1, n);
        end
    end
end
