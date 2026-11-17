struct Trunc
    coeffs::Array{Interval, 1}
    Trunc(a::Array{Interval, 1}) = length(a) == 1+N ? new(a) : error("wrong #coeffs")
end

abs(t::Trunc) = sum(abs(a) for a in t.coeffs)

(^)(t::Trunc, p::Integer) = russian(t, p)

zero(t::Trunc) = Trunc([Interval(0) for k in 0:N])
one(t::Trunc) = Trunc([k == 0 ? Interval(1) : Interval(0) for k in 0:N])

# Vector space essentials
(+)(t::Trunc, u::Trunc) = Trunc(t.coeffs + u.coeffs)
(-)(t::Trunc) = Trunc(-t.coeffs)
(*)(t::Trunc, a::Union{RR, Interval}) = Trunc(t.coeffs * a)

# Vector space derived
(-)(t::Trunc, u::Trunc) = t+(-u)
(/)(t::Trunc, a::Interval) = t*(1/a)
(*)(a::Union{RR, Interval}, t::Trunc) = t*a

# I think that (*)(Array{T, 1}, Number) is defined already by default.
# But since Interval is not a subtype of Number,
# we must define (*)(Array{Interval, 1}, Interval).
(*)(a::Array{Interval, 1}, b::Interval) = [ak*b for ak in a]
(*)(b::Interval, a::Array{Interval, 1}) = a*b

# Convenience
function +(t::Trunc, a::Interval)
    coeffs = copy(t.coeffs)
    coeffs[1] += a
    Trunc(coeffs)
end

function -(t::Trunc, a::Interval)
    coeffs = copy(t.coeffs)
    coeffs[1] -= a
    Trunc(coeffs)
end

(+)(a::Interval, t::Trunc) = t+a
(-)(a::Interval, t::Trunc) = (-t)+a

# Algebra
function *(t::Trunc, u::Trunc)
    coeffs = fill(Interval(0), N+1)
    for j in 0:N
        for k in 0:N-j
            coeffs[1+j+k] += t.coeffs[1+j]*u.coeffs[1+k]
        end
    end
    Trunc(coeffs)
end

# Evaluation and composition of polynomial. Compact, but could be made more efficient?
#(t::Trunc)(x::Union{Interval, Rectangle, Trunc}) = sum(a*x^(k-1) for (k, a) in enumerate(t.coeffs))
function (t::Trunc)(x::Union{Interval, Rectangle, Trunc})
    total = t.coeffs[N+1]*x + t.coeffs[N]
    for k in N-1:-1:1
        total = total*x + t.coeffs[k]
    end
    total
end

function bound_high_order_product(t::Trunc, u::Trunc)
    coeffs = Dict{Int64, Interval}()
    for j in 0:N
        for k in 0:N
            if j+k > N
                a::Interval = t.coeffs[1+j]*u.coeffs[1+k]
                if haskey(coeffs, j+k)
                    coeffs[j+k] += a
                else
                    coeffs[j+k] = a
                end
            end
        end
    end
    sum(abs(a) for (k, a) in coeffs)
end

abs_star(t::Trunc) = sum(abs(t.coeffs[1+k]) for k in 1:N)

diff(t::Trunc) = Trunc([k == N ? Interval(0) : (1+k)*t.coeffs[2+k] for k in 0:N])
