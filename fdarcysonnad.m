function f=fdarcysonnad(ks,D,Q,v)
% function f=fdarcysonnad(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow ---> Colebrook-White 
% 
% EXPLICIT FORMULATION, Sonnad and Goudar (2006)
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
  G = 0.124*Re*Fr + log(0.4587*Re);
  C1 = 0.4587*Re;
  C2 = (G-0.31)^(G/(G+1));
  A = 0.8686*log(C1/C2);
  f = 1/(A*A);
end