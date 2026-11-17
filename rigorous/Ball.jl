struct Ball
    std::StdBall
    dom::Domain
end

# negative
(-)(f::Ball) = Ball(-f.std, f.dom)

# scalar multiply
(*)(f::Ball, i::Interval) = Ball(f.std*i, f.dom)
(*)(i::Interval, f::Ball) = f*i
(/)(f::Ball, i::Interval) = f*(1/i)

(*)(f::Ball, a::Number) = f*Interval(a)
(*)(a::Number, f::Ball) = f*a

# ball add
(+)(f::Ball, g::Ball) = f.dom == g.dom ? Ball(f.std+g.std, f.dom) : error("domain mismatch")
(-)(f::Ball, g::Ball) = f.dom == g.dom ? Ball(f.std-g.std, f.dom) : error("domain mismatch")

# scalar add
(+)(f::Ball, a::Interval) = Ball(f.std+a, f.dom)
(-)(f::Ball, a::Interval) = Ball(f.std-a, f.dom)

# ball mul
(*)(f::Ball, g::Ball) = f.dom == g.dom ? Ball(f.std*g.std, f.dom) : error("domain mismatich")

# ball pow
(^)(f::Ball, n::Integer) = Ball(f.std^n, f.dom)

# ball eval/compose
(f::Ball)(i) = f.std((i - f.dom.c)/f.dom.r)
(f::Ball)(g::Ball) = Ball(f.std((g.std - f.dom.c)/f.dom.r), g.dom)

# diff_compose
diff_compose(f::Ball, g::Ball) = Ball(diff_compose(f.std, (g.std - f.dom.c)/f.dom.r)/f.dom.r, g.dom)

# allow differentiation for composition
function diff(f::Ball)
    df(g::Ball) = diff_compose(f, g)
end

# identity with respect to domain
function identity(f::Ball)
    coeffs = fill(Interval(0), N+1)
    coeffs[1] = f.dom.c
    coeffs[2] = f.dom.r
    Ball(StdBall(Trunc(coeffs), Interval(0), Interval(0)), f.dom)
end

basis_element(f::Ball, k::Integer) = Ball(basis_element(f.std, k), f.dom)

function Ball(p::Poly, h::Number, g::Number)
    c, r = p.cr
    dom = Domain(Interval(c), Interval(r))
    t = Trunc(map(Interval, p.coeffs))
    s = StdBall(t, h, g)
    b = Ball(s, dom)
end

abs(b::Ball) = abs(b.std)

function (m::Array{Interval, 2})(b::Ball, coeff_HH::Interval)
    PP = m*b.std.P.coeffs
    # assume all of G is in the polynomial part
    PG = [hull(m[1+k,:])*balanced(b.std.G) for k in 0:N]
    P = Trunc(PP + PG)
    H = (b.std.H + b.std.G)*abs(coeff_HH)
    G = Interval(0)
    Ball(StdBall(P, H, G), b.dom)
end
