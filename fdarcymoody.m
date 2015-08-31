function f=fdarcymoody(ks,D,Q,v)
% function f=fdarcymoody(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% EXPLICIT FORMULATION, MOODY (1947)
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
  Fr = 2e4*(ks/D);  % Pipe Roughness Factor
  Ft = 1e6/Re; % Flow Turbulence Factor
  f  = 0.0055*(1+(Fr+Ft)^(1/3));
  
end