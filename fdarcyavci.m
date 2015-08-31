function f=fdarcyavci(ks,D,Q,v)
% function f=fdarcyavci(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow ---> Colebrook-White 
% 
% EXPLICIT FORMULATION, Avci & Karagoz (2009)
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
  Fr = ks/D;  % Pipe Roughness Factor
  C1 = 6.4; 
  C2 = (log(Re) - log(1+0.01*Re*Fr*(1+10*(Fr)^0.5)))^2.4;
  f  = C1/C2; 
end
