function f=fdarcymanadilli(ks,D,Q,v)
% function f=fdarcymanadilli(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% EXPLICIT FORMULATION, MANADILLI (1997)
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
  C1 = +95.00/Re^0.983;
  C2 = -96.82/Re;
  Ft = C1+C2;  % Flow Turbulence Factor
  Fr = e/3.70; % Pipe Roughness Factor
  A = log10(Fr+Ft); % this formula is with l(x)
  B = 0.5^2;
  f = B/(A^2);
end