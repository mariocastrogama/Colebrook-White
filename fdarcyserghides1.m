function f=fdarcyserghides1(ks,D,Q,v)
% function f=fdarcyserghides(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% SEMI-EXPLICIT FORMULATION, SERGHIDES 1st formulation (1984)
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
  S1 = -2*log10(e/3.7+12/Re);
  S2 = -2*log10(e/3.7+2.51/Re*S1);
  S3 = -2*log10(e/3.7+2.51/Re*S2);
  A = S1-(S2-S1)^2/(S3-2*S2+S1); 
  B = 1;
  f = B/(A^2);
end