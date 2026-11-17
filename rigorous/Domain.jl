struct Domain
    c::Interval
    r::Interval
    Domain(c::Interval, r::Interval) = r.lo <= 0 ? error("domain radius") : new(c, r)
end

in(z::Union{BigFloat, Interval}, d::Domain) = abs(z-d.c) < d.r
in(z::Rectangle, d::Domain) = abs(z-Rectangle(d.c, 0)) < d.r
(==)(d::Domain, e::Domain) = d.c==e.c && d.r==e.r
