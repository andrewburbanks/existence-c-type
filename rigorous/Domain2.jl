struct Domain2
    c::Tuple{Interval, Interval}
    r::Tuple{Interval, Interval}
    Domain2(c::Tuple{Interval, Interval}, r::Tuple{Interval, Interval}) = r[1].lo > 0 && r[2].lo > 0 ? new(c, r) : error("all radii must be positive")
end

function in(x::Union{BigFloat, Interval}, d::Domain2) 
    abs(x - d.c[1]) < d.r[1] || abs(x - d.c[2]) < d.r[2] 
end

(==)(d1::Domain2, d2::Domain2) = d1.c[1] == d2.c[1] && d1.c[2] == d2.c[2] && d1.r[1] == d2.r[1] && d1.r[2] == d2.r[2]

ψ(d::Domain2, z) = (z .- d.c) ./ d.r
