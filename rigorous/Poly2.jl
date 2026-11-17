struct Poly2
    coeffs::Array{RR, 2}
    Poly2(c) = new([1+j <= size(c)[1] && 1+k <= size(c)[2] && j+k <= N ? c[1+j,1+k] : RR(0) for j=0:N,k=0:N])
end

# unary operators
(-)(p::Poly2) = Poly2(.-p.coeffs)
(+)(p::Poly2) = Poly2(p.coeffs)

# polynomial addition
(+)(p::Poly2, q::Poly2) = Poly2(p.coeffs .+ q.coeffs)
(-)(p::Poly2, q::Poly2) = p + (-q)

# polynomial/number addition
(+)(p::Poly2, q::RR) = p + Poly2(fill(q, (1,1)))
(+)(q::RR, p::Poly2) = p + q
(-)(p::Poly2, q::RR) = p+(-q)
(-)(q::RR, p::Poly2) = -p+q

# scalar multiply
(*)(p::Poly2, a::RR) = Poly2(p.coeffs .* a)
(*)(a::RR, p::Poly2) = p*a 
(/)(p::Poly2, a::RR) = Poly2(p.coeffs ./ a)

# norm
abs(p::Poly2) = sum(abs.(p.coeffs))
abs_star(p::Poly2) = sum(abs.(p.coeffs[2:end]))

# the following might be needed
#absr(p::Poly2) = abs(p)
#absr_star(p::Poly2) = abs_star(p)

# polynomial multiplication
function (*)(p::Poly2, q::Poly2)
    c = fill(RR(0), (N+1,N+1))
    for jp=0:N
        for jq=0:N-jp
            for kp=0:N-jp-jq #triangular coeff
                for kq=0:N-jp-jq-kp #triangular coeff
                    c[1+jp+jq,1+kp+kq] += p.coeffs[1+jp,1+kp]*q.coeffs[1+jq,1+kq]
                end
            end
        end
    end
    Poly2(c)
end

# power of polynomial
(^)(t::Poly2, p::Integer) = (p >= 0) ? russian(t, p) : error("not supported")

# evaluation and composition
function (p::Poly2)(x, y)
    s = x*y*RR(0)  
    for j=N:-1:0
        t = x*y*RR(0)
        for k=N-j:-1:0 #triangular coeff
            t *= y
            t += p.coeffs[1+j,1+k]
        end
        s *= x
        s += t
    end
    s
end

# differentiation 
d1(p::Poly2) = Poly2([(1+j)*p.coeffs[2+j,1+k] for j=0:N-1,k=0:N])
d2(p::Poly2) = Poly2([(1+k)*p.coeffs[1+j,2+k] for j=0:N,k=0:N-1])

# unit vectors given number of elements and index
unit_vector(k::Int64, n::Int64) = [j == k ? RR(1) : RR(0) for j=0:n-1]

# coordinate functions
coordinate_maps(p::Tuple{Poly2, Poly2}) = (Poly2([0 0 0; 1 0 0; 0 0 0]), Poly2([0 1 0; 0 0 0; 0 0 0]))

# coefficient maps; sum of homogeneous polynomials / triangular basis
tri_basis_size_1 = div((N+1)*(N+2), 2) # single map of 2 variables
tri_basis_size_2 = 2*tri_basis_size_1 # pair of maps of 2 variables
tri_basis_size_4 = 4*tri_basis_size_1 # pair of pairs of maps

function tri_phi1(t::Poly2)
    v = fill(RR(0), tri_basis_size_1)
    i = 0
    for j=0:N
        for k=0:N-j
            v[1+i] = t.coeffs[1+j,1+k]
            i += 1
        end
    end
    v
end

function tri_phi1_inv(v::Array{RR, 1})
    c = fill(RR(0), (1+N, 1+N))
    i = 0
    for j=0:N
        for k=0:N-j
            c[1+j,1+k] = v[1+i]
            i += 1
        end
    end
    Poly2(c)
end

tri_phi2(p::Tuple{Poly2, Poly2}) = [tri_phi1(p[1]); tri_phi1(p[2])] 
tri_phi2_inv(c::Vector{RR}) = (tri_phi1_inv(c[1:div(length(c), 2)]), tri_phi1_inv(c[div(length(c), 2)+1:end]))

tri_phi4(p::Tuple{Poly2, Poly2, Poly2, Poly2}) = [tri_phi1(p[1]); tri_phi1(p[2]); tri_phi1(p[3]); tri_phi1(p[4])] 
tri_phi4_inv(c::Vector{RR}) = (tri_phi1_inv(c[1:div(length(c), 4)]), tri_phi1_inv(c[div(length(c), 4)+1:div(length(c), 2)]), tri_phi1_inv(c[div(length(c), 2)+1:3*div(length(c), 4)]), tri_phi1_inv(c[3*div(length(c), 4)+1:end]))

# convenience functions for pairs of polys
(-)(p::Tuple{Poly2, Poly2}, q::Tuple{Poly2, Poly2}) = p .- q
abs(t::Tuple{Poly2, Poly2}) = sum(abs.(t))
