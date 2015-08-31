function f=fdarcyswameerathie(ks,D,Q,v)
% function f=fdarcyswameerathie(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
% Laminar flow   ---> Poiseuille
% Turbulent flow ---> Colebrook-White 
%
% FORMULATION SWAMEE and RATHIE (2007)
% 
% fr: friction factor [n/a] one value at a time not dependent of Reynolds
%     number.
% ks: roughness (average of pipe - channel) [m]
% D : diameter of pipe [m]
% Q : Discharge [m3/s]
% v : cinematic viscosity [m2/s] typical value water 1e-6
% 
% By Mario CASTRO GAMA
% MSc Hydroinformatics
% 2013.01.03
% 
% Requires numre.m for Reynolds number calculation
%
  Re = numre(Q,D,v);
  fr = fdarcyrough(ks,D);
  iRfr = fr^0.5;
  iRfr = 1 / iRfr;
  
  A = 8.0676*D/(ks*Re);
  A2 = A*A;
  temp = 1 - A + A2*(1+0.5756*iRfr);
  temp = temp - A*A2*(1+0.3123*iRfr)*(1+1.4144*iRfr);
  temp = temp + A2*A2*(1+0.2232*iRfr)*(1+0.6663*iRfr)*(1+2.5639*iRfr);
  temp = temp * iRfr;
  f = 1 / (temp * temp);
end