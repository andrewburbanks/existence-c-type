import Base: <, <=, ==, >, >=
import Base: zero, one
import Base: abs
import Base: +, -, *, /, ^

# This struct represents a Poly expanded with respect
# to a disc domain with centre and radius cr = (c, r).

struct Poly
    coeffs::Array{RR, 1}
    cr::Tuple{RR, RR}
    Poly(a::Array{RR, 1}, cr) = length(a) == 1+N ? new(a, cr) : error("wrong #coeffs")
end

abs(t::Poly) = sum(abs(a) for a in t.coeffs)
zero(p::Poly) = Poly(fill(RR(0), N+1), p.cr)
one(p::Poly) = Poly([k==0 ? RR(1) : RR(0) for k in 0:N], p.cr)

(-)(p::Poly) = Poly(-p.coeffs, p.cr)

(*)(p::Poly, a::RR) = Poly(p.coeffs*a, p.cr)
(*)(a::RR, p::Poly) = p*a

(/)(p::Poly, a::RR) = p*(1/a)

function +(p::Poly, a::RR)
    coeffs = copy(p.coeffs)
    coeffs[1] += a
    Poly(coeffs, p.cr)
end

(+)(a::RR, p::Poly) = p+a
(-)(p::Poly, a::RR) = p+(-a)

(+)(p::Poly, q::Poly) = p.cr == q.cr ? Poly(p.coeffs+q.coeffs, p.cr) : error("!")
(-)(p::Poly, q::Poly) = p.cr == q.cr ? Poly(p.coeffs-q.coeffs, p.cr) : error("!")

(*)(a::Array{RR, 1}, b::RR) = [ak*b for ak in a]
(*)(b::RR, a::Array{RR, 1}) = a*b

(/)(a::Array{RR, 1}, b::RR) = [ak/b for ak in a]

function *(p::Poly, q::Poly)
    if p.cr != q.cr error("!") end
    coeffs = fill(RR(0), N+1)
    for j in 0:N
        for k in 0:N
            if j+k <= N
                coeffs[1+j+k] += p.coeffs[1+j]*q.coeffs[1+k]
            end
        end
    end
    Poly(coeffs, p.cr)
end

function identity(p::Poly)
    c, r = p.cr
    coeffs = fill(RR(0), N+1)
    coeffs[1] = c
    coeffs[2] = r
    Poly(coeffs, p.cr)
end

function diff(p::Poly)
    c, r = p.cr
    Poly([k == N+1 ? RR(0) : k*p.coeffs[1+k] for k in 1:N+1], p.cr)/r
end

(^)(t::Poly, p::Integer) = russian(t, p)

# Evaluation and composition of polynomial.
function (t::Poly)(x::Union{RR, Poly, Complex{RR}})
    c, r = t.cr
    u = (x-c)/r
    total = t.coeffs[N+1]*u + t.coeffs[N]
    for k in N-1:-1:1
        total = total*u + t.coeffs[k]
    end
    total
end
