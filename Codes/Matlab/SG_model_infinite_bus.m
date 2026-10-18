% SG_model_infinite_bus.m
% Model of Synchronous Generator withInfinite Bus in Point of Connection
  function dx = SG_model_infinite_bus(t,x,freq,vf0,Pm,V,theta,xd,x1d,x2d,T1d0,T2d0,xq,x1q,x2q,T1q0,T2q0,H,D,n)
% Authors: Alex Stankovic and Alex Saric
% Edits: Explicitly rewriting in terms of id, iq, Pg to avoid confusion.
% -mkt 12/15/2020

   omega0=1;
   omegab=2*pi*freq;
   vd = V*sin(x(1) - theta);
   vq = V*cos(x(1) - theta);   
   if n==6,   
       id = ( x(5) - vq)/x2d;
       iq = (vd - x(6) )/x2q;
       Pg = vd*id + vq*iq;
       dx = [
             omegab*(x(2)-omega0);
             1/H*(Pm - Pg -D*(x(2)-omega0));
             1/T1d0*(-x(3)-(xd-x1d)*id+vf0);
             1/T1q0*(-x(4)+(xq-x1q)*iq);
             1/T2d0*(-x(5)+x(3)-(x1d-x2d)*id);
             1/T2q0*(-x(6)+x(4)+(x1q-x2q)*iq)
            ];
   elseif n==4,
       id = ( x(3) - vq)/x1d;
       iq = (vd - x(4) )/x1q;
       Pg = vd*id + vq*iq;
       dx = [
             omegab*(x(2)-omega0);
             1/H*(Pm-Pg -D*(x(2)-omega0));
             1/T1d0*(-x(3)-(xd-x1d)*id+vf0);
             1/T1q0*(-x(4)+(xq-x1q)*iq)
             ];
   elseif n==3,
       id = ( x(3) - vq)/x1d;
       iq = vd/xq;
       Pg = vd*id + vq*iq;
       dx = [
             omegab*(x(2)-omega0);
             1/H*(Pm-Pg -D*(x(2)-omega0));
             1/T1d0*(-x(3)-(xd-x1d)*id+vf0)
             ];             
   elseif n==2,
       eqp=V+x1d*Y(4,1)/V;
       dx = [
             omegab*(x(2)-omega0);
             1/H*(Pm-V*sin(x(1)-theta)*(eqp-V*cos(x(1)-theta))/x1d-V*cos(x(1)-theta)*(V*sin(x(1)-theta))/x1q)-D*(x(2)-omega0);
             ];
   end
end


