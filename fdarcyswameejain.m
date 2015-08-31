function f=fdarcyswameejain(ks,D,Q,v)
% function f=fdarcyswameejain(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
% based on the formulation of Swamee P.K. & Jain A.K. (1976), 
% Explicit equations for pipe-flow problems,
% Journal of the Hydraulic division, vol 102, No 5, pp 657-664.
%
% Turbulent flow ---> Colebrook-White
% 
% EXPLICIT FORMULATION, SWAMEE AND JAIN (1976)
%
% f friction factor [n/a] one value at a time
% ks roughness (average of pipe - channel) [m]
% D diameter of pipe [m]
% Q Discharge [m3/s]
% v cinematic viscosity [m2/s] typical value water 1e-6
% 
% By Mario CASTRO GAMA
% MSc Hydroinformatics
% 2013.01.03
% 
% Requires numre.m for Reynolds number calculation
%
  Re=numre(Q,D,v); 
  Fr=(ks/D)/3.7;  % Pipe Roughness Factor
  Ft=5.74/(Re^0.9); % Flow Turbulence Factor
  A=log10(Fr+Ft);
  f=0.25/(A^2);
end