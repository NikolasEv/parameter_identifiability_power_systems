% perturb_param.m
% Perturbation of parameters
  function [xdp,x1dp,x2dp,T1d0p,T2d0p,xqp,x1qp,x2qp,T1q0p,T2q0p,Hp,Dp] = perturb_param(interval,xd,x1d,x2d,T1d0,T2d0,xq,x1q,x2q,T1q0,T2q0,H,D)
% Authors: Alex Stankovic and Alex Saric
% Edits: Re-parameterizing to ensure constraints are satisfied after random
% perturbations

  deltax1 = (xd - xq);
  deltax2 = (xq - x1q);
  deltax3 = (x1q - x1d);
  deltax4 = (x1d - x2q);
  deltax5 = (x2q - x2d);
  
  deltaTd = T1d0 - T2d0;
  deltaTq = T1q0 - T2q0;
  
  ponder=(rand(12,1)-0.5)*interval/100;
  
  x2dp = x2d*(1 + ponder(3,1));
  x2qp = x2dp + deltax5*(1 + ponder(8,1));
  x1dp = x2qp + deltax4*(1 + ponder(2,1));
  x1qp = x1dp + deltax3*(1 + ponder(7,1));
  xqp = x1qp + deltax2*(1 + ponder(6,1));
  xdp = xqp + deltax1*(1 + ponder(1,1));
  
  T2d0p = T2d0*(1 + ponder(5,1));
  T1d0p = T2d0p + deltaTd*(1 + ponder(4,1));
  
  T2q0p = T2q0*(1 + ponder(10,1));
  T1q0p = T2q0p + deltaTq*(1 + ponder(9,1));

%   xdp=xd*(1+ponder(1,1));
%   x1dp=x1d*(1+ponder(2,1));
%   x2dp=x2d*(1+ponder(3,1));
%   T1d0p=T1d0*(1+ponder(4,1));
%   T2d0p=T2d0*(1+ponder(5,1));
%   xqp=xq*(1+ponder(6,1));
%   x1qp=x1q*(1+ponder(7,1));
%   x2qp=x2q*(1+ponder(8,1));
%   T1q0p=T1q0*(1+ponder(9,1));
%   T2q0p=T2q0*(1+ponder(10,1));
  Hp=H*(1+ponder(11,1));
  Dp=D*(1+ponder(12,1));

