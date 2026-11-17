struct Trunc2
    coeffs::Array{Interval, 2}
    function Trunc2(a::Array{Interval, 2})
        new([1+j<=size(a)[1] && 1+k<=size(a)[2] && j+k<=N ? a[1+j,1+k] : Interval(0) for j=0:N, k=0:N])
    end
end

# useful constants
zero(t::Trunc2) = Trunc2([Interval(0) for j in 0:N, k in 0:N])
one(t::Trunc2) = Trunc2([j == 0 && k == 0 ? Interval(1) : Interval(0) for j in 0:N, k in 0:N])

# unary operators
(-)(t::Trunc2) = Trunc2(-t.coeffs)
(+)(t::Trunc2) = Trunc2(t.coeffs)

# trunc addition
(+)(t1::Trunc2, t2::Trunc2) = Trunc2(t1.coeffs + t2.coeffs)
(-)(t1::Trunc2, t2::Trunc2) = t1 + (-t2)

# trunc/number addition
function +(t::Trunc2, i::Interval)
    coeffs = copy(t.coeffs)
    coeffs[1] += i
    Trunc2(coeffs)
end
(+)(i::Interval, t::Trunc2) = t + i
(-)(t::Trunc2, i::Interval) = t + (-i)
(-)(i::Interval, t::Trunc2) = (-t) + i

# scalar multiply
# the following line assumes we have defined (*)(::Matrix{Interval}, ::Interval)
(*)(t::Trunc2, a::Union{RR, Interval}) = Trunc2(t.coeffs * a)
(*)(a::Union{RR, Interval}, t::Trunc2) = t * a
(/)(t::Trunc2, a::Union{RR, Interval}) = t * (1/a)

# norm
absr(t::Trunc2) = sum(absr.(t.coeffs))
absr_star(t::Trunc2) = sum(absr.(t.coeffs[2:end]))

# Check this.
abs(t::Trunc2) = sum(abs.(t.coeffs))
abs_star(t::Trunc2) = sum(abs.(t.coeffs[2:end]))

# trunc multiply
function *(t1::Trunc2, t2::Trunc2)
    coeffs = fill(Interval(0, 0), (N+1, N+1))
    for j1 in 0:N
        for j2 in 0:N-j1
            for k1 in 0:N-j1-j2
                for k2 in 0:N-j1-j2-k1
                    coeffs[1+j1+j2, 1+k1+k2] += t1.coeffs[1+j1, 1+k1] * t2.coeffs[1+j2, 1+k2]
                end
            end
        end
    end
    Trunc2(coeffs)
end

# trunc multiple - high-order part
function bound_on_higher_order_multiplication_trunc2(t1::Trunc2, t2::Trunc2)
    c = fill(Interval(0.0), (2N + 1, 2N+1))
    for j1 = 0:N
        for j2 = 0:N
            for k1 = 0:N
                for k2 = 0:N
                    if j1 + j2 + k1 + k2 > N
                        c[1+j1+j2, 1+k1+k2] += t1.coeffs[1+j1, 1+k1] * t2.coeffs[1+j2, 1+k2]
                    end
                end
            end
        end
    end
    #sum(abs.(c))
    sum(absr.(c))
end  

# power
(^)(t::Trunc2, p::Integer) = (p >= 0) ? russian(t, p) : error("not supported")

#evaluation and composition (double horner method)
function (t::Trunc2)(x, y)
    s1 = x*y*Interval(0)  
    for j=N:-1:0
        s2 = x*y*Interval(0)
        for k=N-j:-1:0 #triangular coeff
            s2 *= y
            s2 += t.coeffs[1+j,1+k]
        end
        s1 *= x
        s1 += s2
    end
    s1
end

# differentiation
d1(t::Trunc2) = Trunc2([(1+j)*t.coeffs[2+j,1+k] for j=0:N-1,k=0:N])
d2(t::Trunc2) = Trunc2([(1+k)*t.coeffs[1+j,2+k] for j=0:N,k=0:N-1])

#unit vectors and basis elements
unit_vector(k::Int64, n::Int64) = [j == k ? Interval(1) : Interval(0) for j=0:n-1]

# coefficient maps
function tri_phi1(t::Trunc2)
    v = fill(Interval(0), tri_basis_size_1)
    i = 0
    for j=0:N
        for k=0:N-j
            v[1+i] = t.coeffs[1+j,1+k]
            i += 1
        end
    end
    v
end

function tri_phi1_inv(v::Array{Interval, 1})
    c = fill(Interval(0), (1+N, 1+N))
    i = 0
    for j=0:N
        for k=0:N-j
            c[1+j,1+k] = v[1+i]
            i += 1
        end
    end
    Trunc2(c)
end

tri_phi2(p::Tuple{Trunc2, Trunc2}) = [tri_phi1(p[1]); tri_phi1(p[2])] 
tri_phi2_inv(c::Vector{Interval}) = (tri_phi1_inv(c[1:div(length(c), 2)]), tri_phi1_inv(c[div(length(c), 2)+1:end]))

tri_phi4(p::Tuple{Trunc2, Trunc2, Trunc2, Trunc2}) = [tri_phi1(p[1]); tri_phi1(p[2]); tri_phi1(p[3]); tri_phi1(p[4])] 
tri_phi4_inv(c::Vector{Interval}) = (tri_phi1_inv(c[1:div(length(c), 4)]), tri_phi1_inv(c[div(length(c), 4)+1:div(length(c), 2)]), tri_phi1_inv(c[div(length(c), 2)+1:3*div(length(c), 4)]), tri_phi1_inv(c[3*div(length(c), 4)+1:end]))

#(-)(t1::Tuple{Trunc2, Trunc2}, t2::Tuple{Trunc2, Trunc2}) = (t1[1] - t2[1], t1[2] - t2[2])
