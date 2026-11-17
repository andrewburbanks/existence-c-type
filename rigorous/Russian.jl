# Power
function russian(t, p::Integer)
    if p == 0
        return one(t)
    elseif p == 1
        return t
    end
    
    # ignore trailing zeros in binary expansion
    trailing = trailing_zeros(p) + 1
    p >>= trailing
    
    # repeatedly square
    while (trailing -= 1) > 0
        t *= t
    end
    
    # compute the rest
    u = t
    while p > 0
        trailing = trailing_zeros(p) + 1
        p >>= trailing
        while (trailing -= 1) >= 0
            t *= t
        end
        u *= t
    end
    u
end
