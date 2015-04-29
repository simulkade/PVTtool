% small test on using negflash for bubbledew calculations
clc; clear;
%
vf0 = 0.1; % vapor fraction, zero means bubble point, one means dew point
T0 = 300; %[K]
p0 = 10e5; % [Pa]
fzero(@(p)(vf0-testflashpT(p,T0)), p0)