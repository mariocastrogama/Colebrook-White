function f=fdarcyeck(ks,D,Q,v)
% function f=fdarcyeck(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% EXPLICIT FORMULATION, ECK (1973)
%
% f friction factor [n/a] one value at a time
% ks roughness (average of pipe - channel) [m]
% D diameter of pipe [m]
% Q Discharge [m3/s]
% v cinematic viscosity [m2/s] typical value water 1e-6
% 
% By Mario CASTRO GAMA
% MSc Hydroinformatics
% 2013.12.23
% 
% Requires numre.m for Reynolds number calculation
%
  Re = numre(Q,D,v);
  Ft = 15/Re;         % Pipe Turbulence Factor
  Fr = ks/(3.715*D);  % Pipe Roughness Factor
  A = -2*log10(Fr+Ft);
  f  = 1/(A*A);
  
end