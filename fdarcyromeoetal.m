function f=fdarcyromeoetal(ks,D,Q,v)
% function f=fdarcyromeoetal(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% SEMI-IMPLICIT FORMULATION, ROMEO, ROYO & MONZON (2002)
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
  e = (ks/D);  
  C1 = log10((e/7.7918)^0.9924+(5.3326/(208.815+Re))^0.9345);
  C2 = log10(e/3.827-4.567/Re*C1);
  Ft = -5.0272/Re*C2;  % Flow Turbulence Factor
  Fr = e/3.7065;          % Pipe Roughness Factor
  A = log10(Fr+Ft);
  B = (0.5)^2;
  f = B/(A^2);
end