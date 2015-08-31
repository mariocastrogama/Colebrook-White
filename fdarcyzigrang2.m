function f=fdarcyzigrang2(ks,D,Q,v)
% function f=fdarcyzigrang2(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% SEMI-IMPLICIT FORMULATION, ZIGRANG 2nd formulae (1982)
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
  C1 = 13.0/Re;
  C2 = 5.02/Re;
  Ft = -C2*log10(Fr-C2*log10(Fr+C1));  % Flow Turbulence Factor
  A = log10(Fr+Ft);
  B = (0.5)^2;
  f = B/(A^2);
end