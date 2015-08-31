function f=fdarcyround(ks,D,Q,v)
% function f=fdarcyround(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% EXPLICIT FORMULATION, ROUND (1980)
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
  Fr = 0.135*(ks/D);  % Pipe Roughness Factor
  Ft = 6.5/(Re^1.0); % Flow Turbulence Factor
  A = log10(Fr+Ft);
  B = (-1/1.8)^2;
  f = B/(A^2);
end