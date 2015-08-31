function f=fdarcyhaaland(ks,D,Q,v)
% function f=fdarcyhaaland(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% EXPLICIT FORMULATION HAALAND (1983)
%
% f friction factor [n/a] one value at a time
% ks roughness (average of pipe - channel) [m]
% D diameter of pipe [m]
% Q Discharge [m3/s]
% v cinematic viscosity [m2/s] typical value water 1e-6
% 
% By Mario CASTRO GAMA
% MSc Hydroinformatics
% 2012.12.13
% 
% Requires numre.m for Reynolds number calculation
%
  Re = numre(Q,D,v);
  n  = 1; % 1. liquid, 3. Gas
  Fr = (ks/(3.7*D))^(1.11*n);  % Pipe Roughness Factor
  Ft = (6.9/Re)^n;             % Flow Turbulence Factor
  A = log10(Fr+Ft);
  B = (-1/1.8)^2;
  f = B/(A^2);
end