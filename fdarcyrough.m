function fr=fdarcyrough(ks,D)
% function f=fdarcyrough(ks,D)
% Calculate the friction factor for Darcy-Weisbach equation
% Laminar flow   ---> Poiseuille
% Turbulent flow ---> Colebrook-White 
%
% FORMULATION ONLY FOR ROUGH PIPES Prandtl-Karman as presented by ROUSE H. (1946)
% 
% fr: friction factor [n/a] one value at a time not dependent of Reynolds
%    number.
% ks: roughness (average of pipe - channel) [m]
% D diameter of pipe [m]
% 
% By Mario CASTRO GAMA
% MSc Hydroinformatics
% 2013.01.03
% 
%

  a = (ks/D)/3.70;
  a = log(a);
  fr = 1.325 / (a*a);
end