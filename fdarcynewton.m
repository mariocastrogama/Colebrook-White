function f=fdarcynewton(ks,D,Q,v,tol)
% function f=fdarcynewton(ks,D,Q,v)
% Calculate the friction factor for Darcy-Weisbach equation
% Laminar flow   ---> Poiseuille
% Turbulent flow ---> Colebrook-White 
%
% IMPLICIT FORMULATION - NEWTON RAPHSON [Saldarriaga, 2004]
% 
% f  : friction factor [n/a] one value at a time
% ks : roughness (average of pipe - channel) [m]
% D  : Diameter of pipe [m]
% Q  : Discharge [m3/s]
% v  : cinematic viscosity [m2/s] typical value water 1e-6
% 
% By Mario CASTRO GAMA
% MSc Hydroinformatics
% 2012.12.13
% 
% Requires numre.m for Reynolds number calculation
%
  delta=-100;
  Re=numre(Q,D,v);
  a  = (ks/D)/3.70;
  b  = 2.51/Re;
  if nargin<5;
    tol=1e-9;
  end
  
  if Re<2200; % Laminar Flow
    f = 64/Re;
  else
    f = 0.0001;
    %f = fdarcyswameejain(ks,D,Re); % Explicit approximation by Swamee & Jain (1976)
    while abs(delta) > tol;
      x = 1/f^0.5;
      Fx = -2.0 * log10(a + b*x);
      DFx = (-2/log(10))*(b/(a + b*x));
      delta = (Fx - x)/ (DFx - 1);
      x = x - delta;
      f = 1 / (x*x);
    end
  end
end