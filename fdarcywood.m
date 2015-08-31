function f=fdarcywood(ks,D,Q,v)
% function f=fdarcywood(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% EXPLICIT FORMULATION, WOOD (1966)
%
% f friction factor [n/a] one value at a time
% ks roughness (average of pipe - channel) [m]
% D diameter of pipe [m]
% Q Discharge [m3/s]
% v cinematic viscosity [m2/s] typical value water 1e-6
% 
% By Mario CASTRO GAMA
% MSc Hydroinformatics
% 2013.12.23
% 
% Requires numre.m for Reynolds number calculation
%
  Re = numre(Q,D,v); 
  Fr = ks/D;  % Pipe Roughness Factor
  V  = 1.62*Fr^0.134;
  f  = 0.094*Fr^0.225 + 0.53*Fr + 88*(Fr^0.44)*(Re^(-V));
  
end