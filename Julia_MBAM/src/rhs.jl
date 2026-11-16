function rhs(ps::SG{T}, t, y, dy) where  T<:Real
    H, D, xd, x1d, T1d0, T2d0, xq, x1q, x2q, T1q0, T2q0 = ps()
    vd = V .* sin(y[1] - theta);
    vq = V .* cos(y[1] - theta);   
    id = ( y[5] - vq)/x2d;
    iq = (vd - y[6] )/x2q;
    Pg = vd*id + vq*iq;
    dy[1] = omegab*(y[2]-omega0);
    dy[2] = 1/H*(Pm- Pg -D*(y[2]-omega0));
    dy[3] = 1/T1d0*(-y[3]-(xd-x1d)*id+vf0);
    dy[4] = 1/T1q0*(-y[4]+(xq-x1q)*iq);
    dy[5] = 1/T2d0*(-y[5]+y[3]-(x1d-x2d)*id);
    dy[6] = 1/T2q0*(-y[6]+y[4]+(x1q-x2q)*iq)
    nothing
end

function res(ps::SG{T}, t, y, dy, out) where T<:Real
    H, D, xd, x1d, T1d0, T2d0, xq, x1q, x2q, T1q0, T2q0 = ps()
    vd = V .* sin(y[1] - theta);
    vq = V .* cos(y[1] - theta);   
    id = ( y[5] - vq)/x2d;
    iq = (vd - y[6] )/x2q;
    Pg = vd*id + vq*iq;
    out[1] = omegab*(y[2]-omega0) - dy[1]
    out[2] = (Pm- Pg -D*(y[2]-omega0)) - H * dy[2]
    out[3] = 1/(1 + T1d0)*(-y[3]-(xd-x1d)*id+vf0) - xform(T1d0)*dy[3]
    out[4] = 1/(1 + T1q0)*(-y[4]+(xq-x1q)*iq) - xform(T1q0)*dy[4]
    out[5] = 1/(1 + T2d0)*(-y[5]+y[3]-(x1d-x2d)*id) - xform(T2d0)*dy[5]
    out[6] = 1/(1 + T2q0)*(-y[6]+y[4]+(x1q-x2q)*iq) - xform(T2q0)*dy[6]
    nothing
end
