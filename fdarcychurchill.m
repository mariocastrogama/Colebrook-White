function f=fdarcychurchill(ks,D,Q,v)
% function f=fdarcychurchill(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% EXPLICIT FORMULATION, CHURCHILL (1973)
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
  Fr = (ks/D)/3.7;  % Pipe Roughness Factor
  Ft = 7.0/(Re^0.9); % Flow Turbulence Factor
  A = log10(Fr+Ft);
  B = (1/2)^2;
  f = B/(A^2);
end