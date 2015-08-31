function f=fdarcychen(ks,D,Q,v)
% function f=fdarcychen(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% SEMI-IMPLICIT FORMULATION, CHEN (1979)
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
  e  = ks/D;  
  C1 = 1/2.8257*e^1.1098;
  C2 = 5.8506/Re^0.8961;
  Ft = -5.0452/Re*log10(C1+C2);  % Flow Turbulence Factor
  Fr = e/3.7065; % Pipe Roughness Factor
  A = log10(Fr+Ft);
  B = (0.5)^2;
  f = B/(A^2);
end