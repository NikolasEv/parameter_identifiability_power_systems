function obs(ps::SG{T}, t, y) where T <: Real
    return T[yi for yi in y]
end
