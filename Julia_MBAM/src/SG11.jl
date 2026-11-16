@parameterspace mutable struct SG
    # deltax1 = (xd - xq), exp
    # deltax2 = (xq - x1q), exp
    # deltax3 = (x1q - x1d), exp
    # deltax4 = (x1d - x2q), exp
    # deltax5 = 0, constant
    # deltaTd = T1d0 - T2d0, exp
    # deltaTq = T1q0 - T2q0, exp
    # H = H, exp
    # D = D, exp
    # T2d0 = T2d0, exp
    # T2q0 = T2q0, exp
    # x2d = x2d, exp

    # Reordering the parameters
    # The commented section are the original ordering
    H = H, exp
    D = D, exp
    deltax1 = (xd - xq), exp
    deltax2 = (xq - x1q), exp
    deltax3 = (x1q - x1d), exp
    deltax4 = (x1d - x2q), exp
    x2d = x2d, exp
    deltaTd = T1d0 - T2d0, exp
    deltaTq = T1q0 - T2q0, exp
    T2d0 = T2d0, exp
    T2q0 = T2q0, exp
    deltax5 = 0, constant    
end

function (ps::SG{T})() where T <: Real
    x2d = ps.x2d
    T2q0 = ps.T2q0
    T2d0 = ps.T2d0
    H = ps.H
    D = ps.D
    x2q = x2d + ps.deltax5
    x1d = x2q + ps.deltax4
    x1q = x1d + ps.deltax3
    xq = x1q + ps.deltax2
    xd = xq + ps.deltax1
    T1d0 = T2d0 + ps.deltaTd
    T1q0 = T2q0 + ps.deltaTq

    return T[H, D, xd, x1d, T1d0, T2d0, xq, x1q, x2q, T1q0, T2q0]
end

differential_vars = [true, true, true, true, true, true]
