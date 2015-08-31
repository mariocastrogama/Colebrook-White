function f=fdarcypapaevangelouetal(ks,D,Q,v)
% function f=fdarcypapaevangelouetal(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% EXPLICIT FORMULATION, PAPAEVANGELOU et al. (2010)
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
  numer = 0.2479 - 0.0000947*(7 - log10(Re))^4;
  Ft = 7.366/Re^0.9142;         % Pipe Turbulence Factor
  Fr = ks/(3.615*D);            % Pipe Roughness Factor
  denom = (log10(Fr+Ft))^2;
  f  = numer/denom;  
end