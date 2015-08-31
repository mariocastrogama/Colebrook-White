function f=fdarcybuzzelli(ks,D,Q,v)
% function f=fdarcybuzzelli(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% SEMI-EXPLICIT FORMULATION, BUZZELLI (2008)
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
  e  = (ks/D)/3.7;  
  %B1 = -2*log10(e+(1.9/Re)*(-2*log(e)));
  B1 = 0.7741*log(Re)-1.41; % This is Ln(x)
  B2 = e*Re+2.51*B1;
  A = B1-(B1+2*log10(B2/Re))/(1+2.18/B2); 
  B = 1;
  f = B/(A*A);
end