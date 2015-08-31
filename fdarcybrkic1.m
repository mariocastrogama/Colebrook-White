function f=fdarcybrkic1(ks,D,Q,v)
% function f=fdarcybrkic1(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% EXPLICIT FORMULATION, BRKIC (2011) First Formulation
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
  B  = log(Re/(1.816*log(1.1*Re/(log(1+1.1*Re))))); % log is ln in matlab
  Fr = ks/(3.71*D);            % Pipe Roughness Factor
  Ft = 10^(-0.4343*B);         % Pipe Turbulence Factor
  A = -2*log10(Fr+Ft);
  f  = 1/(A*A);  
end