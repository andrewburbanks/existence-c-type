struct StdBall
    P::Trunc
    H::Interval
    G::Interval
    
    # we enforce the intervals h, g to be of the form [0, b].
    StdBall(p::Trunc, h::RR, g::RR) = new(p, Interval(0, h), Interval(0, g))
    StdBall(p::Trunc, h::Interval, g::Interval) = new(p, Interval(0, h.hi), Interval(0, g.hi))
end

zero(b::StdBall) = StdBall(zero(b.P), Interval(0), Interval(0))
one(b::StdBall) = StdBall(one(b.P), Interval(0), Interval(0))

# Abs
abs(b::StdBall) = nonnegative(abs(b.P) + b.H + balanced(b.G))

# Vector space operations

# Negative
(-)(b::StdBall) = StdBall(-b.P, b.H, b.G)

# Scalar multiples
(*)(b::StdBall, i::Interval) = StdBall(b.P*i, b.H*abs(i), b.G*abs(i))
(/)(b::StdBall, i::Interval) = b*(1/i)

# Add and subtract
(+)(b::StdBall, c::StdBall) = StdBall(b.P+c.P, b.H+c.H, b.G+c.G)
(-)(b::StdBall, c::StdBall) = b+(-c)

# Constant add and subtract
(+)(b::StdBall, i::Interval) = StdBall(b.P+i, b.H, b.G)
(-)(b::StdBall, i::Interval) = StdBall(b.P-i, b.H, b.G)

# Derived
(*)(i::Interval, b::StdBall) = b*i
(+)(i::Interval, b::StdBall) = b+i
(-)(i::Interval, b::StdBall) = (-b)+i

function *(b::StdBall, c::StdBall)
    PP_P = b.P*c.P
    PP_H = bound_high_order_product(b.P, c.P)
    PH = abs(b.P)*c.H
    PG = abs(b.P)*c.G

    HP = b.H*abs(c.P)
    HH = b.H*c.H
    HG = b.H*c.G
    
    GP = b.G*abs(c.P)
    GH = b.G*c.H
    GG = b.G*c.G
    
    P = PP_P
    H = PP_H + PH + HP + HH + GH + HG
    G = PG + GP + GG
    
    StdBall(P, H, G)
end

(^)(t::StdBall, p::Integer) = russian(t, p)

function abs_star(g::StdBall)
    x = abs_star(g.P) + Interval(0, g.G.hi) + Interval(-g.G.hi, g.G.hi)
    Interval(max(0, x.lo), x.hi)
end

# function ball compose

function (f::StdBall)(g::StdBall)
    if ~(f.H == f.G == Interval(0)) && ~(abs(g) in Interval(0, 1))
        error("composition not defined")
    end
    ball_1 = f.P(g)
    ball_3 = StdBall(zero(f.P), Interval(0), f.G)
    g_norm = abs(g)^(N+1)
    g_star_norm = abs_star(g)^(N+1)
    ball_2 = StdBall(zero(f.P), f.H*g_star_norm, f.H*(g_norm - g_star_norm))
    ball_1 + ball_2 + ball_3
end

function (f::StdBall)(g) #could be g::AbstractRigorousScalar
    if ~(f.H == f.G == Interval(0)) && ~(abs(g) in Interval(0, 1))
        error("composition not defined")
    end
    g_norm = abs(g)^(N + 1)
    r = f.H * g_norm + f.G
    f.P(g) + Interval(-r.hi, r.hi)
end

# could absorb this into earlier definition
#(t::Trunc)(x::StdBall) = sum(a*x^(k-1) for (k, a) in enumerate(t.coeffs))
function (t::Trunc)(x::StdBall)
    total = t.coeffs[N+1]*x + t.coeffs[N]
    for k in N-1:-1:1
        total = total*x + t.coeffs[k]
    end
    total
end

function diff_compose(f::StdBall, g::StdBall)
    if (f.H.hi > 0 || f.G.hi > 0) && ~(abs(g) in Interval(0, 1))
        error("diff composition not defined")
    end
    ball_1 = diff(f.P)(g)
    g_norm = abs(g)
    if f.G == f.H == Interval(0)
        norm_sup_H = norm_sup_G = 0
    else
        norm_sup_H = bound_quasi_power(N, g_norm)
        norm_sup_G = bound_quasi_power(0, g_norm)
    end
    # some sacrifice is made here! absorbing into general term!
    # an alternative is to absord into both high-order term and
    # the final polynomial coefficient.
    ball_2 = StdBall(zero(f.P), Interval(0), f.H * norm_sup_H)
    ball_3 = StdBall(zero(f.P), Interval(0), f.G * norm_sup_G)
    ball_1 + ball_2 + ball_3
end

diff_compose(f::StdBall, t::Trunc) = diff_compose(f, StdBall(t, 0, 0))

function basis_element(b::StdBall, k::Integer)
    if k <= N
        return StdBall(Trunc(unit_vector(Interval, k)), 0, 0)
    else
        return StdBall(Trunc(fill(Interval(0), N+1)), 1, 0)
    end
end
