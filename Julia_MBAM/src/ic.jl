# function ic(ps::SG{T}) where T<:Real
#     # Values from Nikos
#     return T[0.5, 0.98, 2.13, 0.02, 1.98, 0.02]

#     # These values from MAT file
#     # return T[1.1764054892711628,
#     #          0.9873199180907988,
#     #          1.2778388802284077,
#     #          0.34302809418086605,
#     #          0.9045839327660419,
#     #          0.8426059099927368]
# end


function ic(ps::SG{T}) where T<:Real
    y = T[0.5, 0.98, 2.13, 0.02, 1.98, 0.02]
    H, D, xd, x1d, T1d0, T2d0, xq, x1q, x2q, T1q0, T2q0 = ps()
    vd = V .* sin(y[1] - theta);
    vq = V .* cos(y[1] - theta);   
    id = ( y[5] - vq)/x2d;
    iq = (vd - y[6] )/x2q;
    Pg = vd*id + vq*iq;
    dy = zeros(T, 6)
    dy[1] = omegab*(y[2]-omega0);
    dy[2] = 1/H*(Pm- Pg -D*(y[2]-omega0));
    dy[3] = 1/T1d0*(-y[3]-(xd-x1d)*id+vf0);
    dy[4] = 1/T1q0*(-y[4]+(xq-x1q)*iq);
    dy[5] = 1/T2d0*(-y[5]+y[3]-(x1d-x2d)*id);
    dy[6] = 1/T2q0*(-y[6]+y[4]+(x1q-x2q)*iq)

    return y, dy
end
