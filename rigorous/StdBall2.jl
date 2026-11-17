struct StdBall2
    t::Trunc2
    H::Interval
    E::Interval
    StdBall2(t::Trunc2, h::RR, e::RR) = (h >= 0 && e >= 0 ? new(t, Interval(0, h), Interval(0, e)) : error("negative norm"))
    StdBall2(t::Trunc2, h::Interval, e::Interval) = (h.hi >= 0 && e.hi >= 0 ? new(t, Interval(0, h.hi), Interval(0, e.hi)) : error("negative norm"))
end

zero(s::StdBall2) = StdBall2(zero(s.t), Interval(0), Interval(0))
one(s::StdBall2) = StdBall2(one(s.t), Interval(0), Interval(0))

(+)(s1::StdBall2, s2::StdBall2) = StdBall2(s1.t + s2.t, s1.H + s2.H, s1.E + s2.E)
(+)(s::StdBall2, i::Interval) = StdBall2(s.t + i, s.H, s.E)
(+)(i::Interval, s::StdBall2)  = s + i

(-)(s::StdBall2) = StdBall2(-s.t, s.H, s.E)
(-)(s1::StdBall2, s2::StdBall2) = s1 + (-s2)
(-)(s::StdBall2, i::Interval) = StdBall2(s.t - i, s.H, s.E)
(-)(i::Interval, s::StdBall2)  = (-s) + i

(*)(s::StdBall2, i::Interval) = StdBall2(s.t * i, s.H * abs(i), s.E * abs(i))
(*)(i::Interval, s::StdBall2)  = s * i
(*)(s::StdBall2, a::RR) = s*Interval(a)
(*)(a::RR, s::StdBall2) = s*a

(/)(s::StdBall2, i::Interval) = s * (1/i)
(/)(s::StdBall2, x::RR) = s * (1/x)

function absr(s::StdBall2)
    x = absr(s.t) + s.H + balanced(s.E)
    Interval(max(0, x.lo), x.hi)
end
function absr_star(s::StdBall2)
    x = absr_star(s.t) + Interval(0, s.H.hi) + balanced(s.E.hi)
    Interval(max(0, x.lo), x.hi)
end

# Check this. Ensure that Trunc2 provides abs.
#ADB ---> abs(s::StdBall2) = absr(s)
function abs(s::StdBall2)
    x = abs(s.t) + s.H + balanced(s.E)
    Interval(max(0, x.lo), x.hi)
end
function abs_star(s::StdBall2)
    x = abs_star(s.t) + Interval(0, s.H.hi) + balanced(s.E.hi)
    Interval(max(0, x.lo), x.hi)
end

function *(f::StdBall2, g::StdBall2)
    PP_P = f.t * g.t
    
    PP_h = bound_on_higher_order_multiplication_trunc2(f.t, g.t)
    PH_h = absr(f.t) * g.H
    HP_h = f.H * absr(g.t)
    HH_h = f.H * g.H
    EH_h = f.E * g.H
    HE_h = f.H * g.E
    h = PP_h + PH_h + HP_h + HH_h + EH_h + HE_h

    PE_e = absr(f.t) * g.E
    EP_e = f.E * absr(g.t)
    EE_e = f.E * g.E
    e = PE_e + EP_e + EE_e
    
    StdBall2(PP_P, h, e)
end

(^)(s::StdBall2, p::Integer) = russian(s, p)

# evaluation at intervals
function (f::StdBall2)(x::Interval, y::Interval)
    @assert abs(x).hi < 1 && abs(y).hi < 1  
    f.t(x, y) + balanced(f.H * hull(abs(x), abs(y))^(N+1) + f.E)
end
# we may need evaluation at rectangles, using balanced rectangle?
function (f::StdBall2)(x::Rectangle, y::Rectangle)
    @assert abs(x).hi < 1 && abs(y).hi < 1
    # the following line needs to produce a single rectangle
    # that contains all possible values f(x,y)=f_P(x,y)+f_H(x,y)+f_E(x,y)
    # where we know that |x|,|y|<1, we have f_P and we know that
    # ||f_H\\<=v_H and f_H is high-order (terms j+k>N) and ||f_E||<=v_E.
    # Note: we need to make a rectangle that covers the disc of all
    # possible values. We do this by making a balanced interval and
    # using it to create a balanced square of the form [-a,+a]x[-a,+a].
    I = balanced(f.H * hull(abs(x), abs(y))^(N+1) + f.E)
    f.t(x, y) + Rectangle(I, I)
end
    
function composition(f::StdBall2, g::StdBall2, h::StdBall2)
    
    @assert absr(g).hi <= 1 && absr(h).hi <= 1
    
    # fP(g, h):
    ball1 = f.t(g,h)
    
    # fH(g, h):
    h_bound = f.H * hull(absr_star(g), absr_star(h))^(N+1)
    bound = Interval(0.0)
    for i = 0:N+1 
        j = N+1-i
        bound = hull(bound, (absr(g)^i * absr(h)^j - absr_star(g)^i * absr_star(h)^j))
    end
    e_bound = f.H * bound
    ball2 = StdBall2(zero(f.t), h_bound, e_bound) 
    
    # fE(g, h):
    ball3 = StdBall2(zero(f.t), Interval(0.0), f.E) 
    
    BALL = ball1 + ball2 + ball3
    BALL
end
(f::StdBall2)(g::StdBall2, h::StdBall2) = composition(f, g, h)

function d1composition(f::StdBall2, g::StdBall2, h::StdBall2)
    c = absr(g)
    d = absr(h)
    e0_bound = bound_quasi_power_2(c, d, 0)
    eN_bound = bound_quasi_power_2(c, d, N)

    # fP'(g, h):
    ball1 = d1(f.t)(g,h)
    # fH'(g, h):
    ball2 = StdBall2(zero(f.t), Interval(0), f.H*eN_bound)
    # fE'(g, h):
    ball3 = StdBall2(zero(f.t), Interval(0), f.E*e0_bound)
    ball1 + ball2 + ball3
end
        
function d2composition(f::StdBall2, g::StdBall2, h::StdBall2)
    c = absr(g)
    d = absr(h)
    e0_bound = bound_quasi_power_2(d, c, 0)
    eN_bound = bound_quasi_power_2(d, c, N)

    ball1 = d2(f.t)(g,h)
    ball2 = StdBall2(zero(f.t), Interval(0), f.H*eN_bound)
    ball3 = StdBall2(zero(f.t), Interval(0), f.E*e0_bound)
    ball1 + ball2 + ball3
end

function d1composition(f::StdBall2, x::Interval, y::Interval)
    c = absr(x)
    d = absr(y)
    e0_bound = bound_quasi_power_2(c, d, 0)
    eN_bound = bound_quasi_power_2(c, d, N)

    I1 = d1(f.t)(x,y)
    I2 = f.H*eN_bound
    I3 = f.E*e0_bound
    I1 + balanced(I2 + I3)
end
        
function d2composition(f::StdBall2, x::Interval, y::Interval)
    c = absr(x)
    d = absr(y)
    e0_bound = bound_quasi_power_2(d, c, 0)
    eN_bound = bound_quasi_power_2(d, c, N)

    I1 = d2(f.t)(x,y)
    I2 = f.H*eN_bound
    I3 = f.E*e0_bound
    I1 + balanced(I2 + I3)
end
