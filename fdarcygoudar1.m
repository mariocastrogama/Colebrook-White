function f=fdarcygoudar1(ks,D,Q,v)
% function f=fdarcygoudar1(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow ---> Colebrook-White 
% 
% EXPLICIT FORMULATION, Goudar C.T. and J.R. Sonnad (2008) First formulation
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
  a = 2/log(10);
  b = ks/D/3.7;  % Pipe Roughness Factor
  d = log(10)/5.02*Re;
  s = b*d + log(d);
  q = s^(s/(s+1));  
  g = b*d + log(d/q);
  z = log(q/g);
  dLA = (g/(g+1))*z;
  f = a*(log(d/q)+dLA);
  f = 1/(f*f);
end