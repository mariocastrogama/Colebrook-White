function f=fdarcybarr(ks,D,Q,v)
% function f=fdarcybarr(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% EXPLICIT FORMULATION, BARR (1981)
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
  Fr = (ks/D)/3.70;  % Pipe Roughness Factor
  C1 = 4.518*log10(Re/7.0); 
  C2 = Re*(1+(Re^0.52)/29*(ks/D)^0.7);
  Ft=C1/C2;         % Flow Turbulence Factor
  A = log10(Fr+Ft);
  B = (0.5)^2;
  f = B/(A^2);
end