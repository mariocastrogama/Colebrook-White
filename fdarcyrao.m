function f=fdarcyrao(ks,D,Q,v)
% function f=fdarcyrao(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow ---> Colebrook-White 
% 
% EXPLICIT FORMULATION, Rao & Kumar (2007)
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
  
  C1 = 1/(2*Fr);
  A  = 0.444/Re+0.135; %(0.444 + 0.135*Re)/Re
  B  = log(Re/6.5);
  B  = -0.33*B*B;
  PhiRe  = 1 - 0.55*exp(B);
  C2 = A*PhiRe;
  f  = 2*log10(C1/C2);
  f  = 1/(f*f);
end