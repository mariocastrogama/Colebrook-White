function f=fdarcysalmashietal(ks,D,Q,v)
% function f=fdarcyslamashietal(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% EXPLICIT FORMULATION, SALMASHI, KHATIBI AND GHORBANI (2012)
%
% f friction factor [n/a] one value at a time
% ks roughness (average of pipe - channel) [m]
% D diameter of pipe [m]
% Q Discharge [m3/s]
% v cinematic viscosity [m2/s] typical value water 1e-6
% 
% By Mario CASTRO GAMA
% MSc Hydroinformatics
% 2014.01.07
% 
% Requires numre.m for Reynolds number calculation
%
  Re = numre(Q,D,v); 
  Fr = (ks/D);  % Pipe Roughness Factor
  temp1 = -11.764*Fr - log(2*Re);
  temp2 = -2.567 + 9.065 / (Re - Fr);
  temp1 = exp(temp1);
  temp2 = exp(temp2);
  f = -0.0575 + Fr + temp1 + temp2;
end