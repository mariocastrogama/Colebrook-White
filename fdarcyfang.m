function f=fdarcyfang(ks,D,Q,v)
% function f=fdarcyfang(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% SEMI-IMPLICIT FORMULATION, FANG (2011)
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
  C1 = -60.525/Re^1.1105;
  C2 = +56.291/Re^1.0712;
  Ft = C1+C2;  % Flow Turbulence Factor
  Fr = 0.234*e^1.1007; % Pipe Roughness Factor
  A = log(Fr+Ft); % this formula is with ln(x)
  B = 1.613;
  f = B/(A^2);
end