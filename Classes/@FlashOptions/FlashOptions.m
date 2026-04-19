classdef FlashOptions
% FlashOptions  Convergence settings for flash and stability calculations.
%
%   opts = FlashOptions()
%
%   All properties have sensible defaults; override as needed:
%     opts.accuracy    = 1e-9;   % tighter tolerance for VLE flash
%     opts.iteration   = 200;    % more iterations for difficult systems
%
%   Properties used by vleflash / vleflashnegative / lleflash:
%     accuracy    - convergence tolerance on mole fraction and Rachford-Rice
%     iteration   - maximum number of successive substitution iterations
%
%   Properties used by stabilityTest / stabilityLLETest:
%     trivialSolutionMaxError - threshold below which solution is trivial
%     convergenceMaxError     - tolerance for convergence of stability loop
%     maxIteration            - maximum stability test iterations

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
        accuracy                = 1e-7
        iteration               = 100
        trivialSolutionMaxError = 1e-5
        convergenceMaxError     = 1e-10
        maxIteration            = 50
    end

    methods
        function obj = FlashOptions()
        end
    end
end
