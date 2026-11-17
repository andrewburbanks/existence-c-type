# Polynomials representing power series expanded with respect to a polydisc domain

struct Poly2D
    c::Array{RR, 1}
    r::Array{RR, 1}
    p::Poly2
end

same_domain(g::Poly2D, f::Poly2D) = (g.c == f.c && g.r == f.r)

(+)(g::Poly2D, f::Poly2D) = same_domain(g, f) ? Poly2D(g.c, g.r, g.p + f.p) : error("Domain mismatch")
(-)(g::Poly2D, f::Poly2D) = same_domain(g, f) ? Poly2D(g.c, g.r, g.p - f.p) : error("Domain mismatch")

(+)(g::Poly2D, a::RR) = Poly2D(g.c, g.r, g.p + a)
(+)(a::RR, g::Poly2D) = g+a

(-)(g::Poly2D) = Poly2D(g.c, g.r, -g.p)
(-)(g::Poly2D, f::RR) = g+(-f)
(-)(a::RR, g::Poly2D) = -g+a

(*)(g::Poly2D, a::RR) = Poly2D(g.c, g.r, g.p * a)
(*)(a::RR, g::Poly2D) = g*a 

(/)(g::Poly2D, a::RR) = Poly2D(g.c, g.r, g.p / a)

abs(g::Poly2D) = abs(g.p)

#norm for a pair of Poly2's
abs(gf::Tuple{Poly2D, Poly2D}) = sum(abs.(gf))

#subtraction for a pair of pairs of Poly2's
(-)(gf1::Tuple{Poly2D, Poly2D}, gf2::Tuple{Poly2D, Poly2D}) = Tuple{Poly2D, Poly2D}(gf1 .- gf2)

# banach algebra product
(*)(g::Poly2D, f::Poly2D) = same_domain(g, f) ? Poly2D(g.c, g.r, g.p * f.p) : error("Domain mismatch")

# evaluation and composition
(g::Poly2D)(x::Union{RR, Complex{RR}}, y::Union{RR, Complex{RR}}) = g.p((([x, y].-g.c)./g.r)...)
# could provide also a Poly2 version, analogous to above.

# take care with this (f o psi) o (g o psi) = (f o psi o g) o psi
# in fact, this ought to propagate the inner domain into the
# resulting poly2d, using the outer one in the expression.
(h::Poly2D)(g::Poly2D, f::Poly2D) = same_domain(g, f) ? Poly2D(g.c, g.r, h.p( (([g.p, f.p].-h.c) ./h.r)...)) : error("Domain mismatch")

# partial derivatives
d1(g::Poly2D) = Poly2D(g.c, g.r, d1(g.p))/g.r[1]
d2(g::Poly2D) = Poly2D(g.c, g.r, d2(g.p))/g.r[2]

# coordinate functions
function coordinate_maps(gf::Tuple{Poly2D, Poly2D})
    g, f = gf
    @assert(same_domain(g, f))
    (cX, cY), (rX, rY) = g.c, g.r
    X_omega = Poly2D(g.c, g.r, Poly2([cX 0; rX 0]))
    Y_omega = Poly2D(g.c, g.r, Poly2([cY rY; 0 0]))
    [X_omega, Y_omega]
end

tri_phi1omega(g::Poly2D) = tri_phi1(g.p) 
tri_phi1omega_inv(c) = Poly2D(domain..., tri_phi1_inv(c))

#for 2 maps of 2 variables
tri_phi2omega(gf::Tuple{Poly2D, Poly2D}) = [tri_phi1omega(gf[1]); tri_phi1omega(gf[2])] 
tri_phi2omega_inv(c) = (tri_phi1omega_inv(c[1:div(length(c), 2)]), tri_phi1omega_inv(c[div(length(c), 2)+1:end]))
