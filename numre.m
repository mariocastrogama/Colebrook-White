function Re=numre(Q,D,v)
% function Re=numre(Q,D,v)
% Calculate the Adimensional Reynolds number of a given flow in a circular pipe.
%
%  Re: Measures the relation between Inertial forces and viscous forces in a
%  Fluid [adimensional].
%
%  Q: Discharge [m3/s,ft^3/s]
%
%  D: Conduit Diameter of pipe [m,ft]. 
%
%  v: viscosity of the fluid [m^2/s,ft^2/s].
% 
% this approach shows the simplified version of the
% Reynolds number usually used in water distribution systems due to 
% the low compressibility of water insisde a pipe or in a channel.
%
  Re=(4*Q)/(pi*D*v);
end