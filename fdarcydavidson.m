function f=fdarcydavidson(ks,D,Q,v)
% function f=fdarcydavidson(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
%
% Turbulent flow --_> Colebrook-White 
% 
% EXPLICIT FORMULATION, DAVIDSON et al. (1999) 14th Formulation
% 
% f friction factor [n/a] one value at a time
% ks roughness (average of pipe - channel) [m]
% D diameter of pipe [m]
% Q Discharge [m3/s]
% v cinematic viscosity [m2/s] typical value water 1e-6
% 
% By Mario CASTRO GAMA
% MSc Hydroinformatics
% 2014.01.03
% 
% Requires numre.m for Reynolds number calculation
%
% Note: the equation has a limited range of aplicability
%       Re [10^5, 10^6]
%       ks/D [0.001, 0.01]
%
  Re = numre(Q,D,v);  
  x1 = 1000*ks/D;     % Pipe Roughness Factor
  x2 = Re/100000;     % Pipe Turbulence Factor
  coef=[1.222995307e-5, -2.242748136e-10, -2.482162347e-4, +9.286977109e-6, 3.645038671e-2,...
        -1.18044694e-3, -0.3849323423, +6.5984401765e-2, 2.522401137, 6.471827292e-4,...
        -1.776888826e-2, +0.1829816121, -0.9369530943, -0.3698214152];
  vars=[x1^6, x1^5*x2^5, x1^5, x1^3*x2^3, x1^3,...
        x1^2*x2^2, x1^2, x1*x2, x1, x2^4,...
        x2^3, x2^2, x1, 1];
  y = coef*vars';
  f  = (y/10)*(0.0385035-0.0199435)+0.0199435;  
end