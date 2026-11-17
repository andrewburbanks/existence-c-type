struct Ball2
    dom::Domain2
    std::StdBall2
end

same_domain(f::Ball2, g::Ball2) = f.dom.c == g.dom.c && f.dom.r == g.dom.r

zero(f::Ball2) = Ball2(f.dom, StdBall2(Trunc2([0]), 0.0, 0.0))
one(f::Ball2) = Ball2(f.dom, StdBall2(Trunc2([1]), 0.0, 0.0))

(+)(f::Ball2, x::Interval) = Ball2(f.dom, f.std + x)
(+)(x::Interval, f::Ball2) = f + x
(+)(f::Ball2, g::Ball2) = same_domain(f, g) ? Ball2(f.dom, f.std + g.std) : error("domain mismatch")

(-)(f::Ball2) = Ball2(f.dom, - f.std)
(-)(f::Ball2, x::Interval) = Ball2(f.dom, f.std - x)
(-)(f::Ball2, g::Ball2) = f + (-g)

(*)(f::Ball2, x::Interval) = Ball2(f.dom, f.std * x)
(*)(x::Interval, f::Ball2) = f * x
(*)(f::Ball2, x::RR) = f*Interval(x)
(*)(x::RR, f::Ball2) = f*x
(*)(f::Ball2, g::Ball2) = same_domain(f, g) ? Ball2(f.dom, f.std * g.std) : error("domain mismatch")

(/)(f::Ball2, x::Interval) = Ball2(f.dom, f.std / x)

(^)(f::Ball2, n::Integer) = Ball2(f.dom, f.std^n)

# missing abs-type functions

#composition
(f::Ball2)(g::Ball2, h::Ball2) = same_domain(g, h) ? Ball2(g.dom, f.std(ψ(f.dom, (g.std, h.std))...)) : error("bad domains")

#evaluation of a Ball2 with numbers or intervals as inputs
(f::Ball2)(x::Rectangle, y::Rectangle) = f.std(ψ(f.dom, (x, y))...) #testing
(f::Ball2)(x::Interval, y::Interval) = f.std(ψ(f.dom, (x, y))...)
(f::Ball2)(x::RR, y::RR) = f(Interval(x), Interval(y)) #can we avoid needing this?

function d1composition(f::Ball2, g::Ball2, h::Ball2)
    Ball2(g.dom, d1composition(f.std, ψ(f.dom, (g.std, h.std))...)/f.dom.r[1])
end

function d2composition(f::Ball2, g::Ball2, h::Ball2)
    Ball2(g.dom, d2composition(f.std, ψ(f.dom, (g.std, h.std))...)/f.dom.r[2])
end

# the following needs to be replaced with
# an efficient implementation of derivative-eval

function d1composition(f::Ball2, x::Interval, y::Interval)
    d1composition(f.std, ψ(f.dom, (x, y))...)/f.dom.r[1]
end

function d2composition(f::Ball2, x::Interval, y::Interval)
    d2composition(f.std, ψ(f.dom, (x, y))...)/f.dom.r[2]
end

# in the following, the intention was that (g, h) be Ball2s,
# however, we also need derivative evaluation as well as
# derivative composition.

# we should perhaps restrict the types of (g, h) in the
# following to be ::Ball2, and then provide a separate
# version for ::Interval?

d1(f::Ball2) = (g, h) -> d1composition(f, g, h)
d2(f::Ball2) = (g, h) -> d2composition(f, g, h)

abs(f::Ball2) = abs(f.std)
absr(f::Ball2) = absr(f.std)

# We need to adapt the following code from Ball to Ball2:
#function (m::Array{Interval, 2})(t::Tuple{Ball2, Ball2}, coeff_HH::Interval)
#    PP = m*b.std.P.coeffs
#    # assume all of G is in the polynomial part
#    PG = [hull(m[1+k,:])*balanced(b.std.G) for k in 0:N] # takes one row of m and multiply by balanced norm G
#    P = Trunc(PP + PG)
#    # assume all of G is in the high-order part
#    H = (b.std.H + b.std.G)*abs(coeff_HH)
#    G = Interval(0)
#    Ball2(b.dom, StdBall2(P, H, G))
#end
