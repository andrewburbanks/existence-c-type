struct Interval
    lo::RR
    hi::RR
    Interval(lo, hi) = lo > hi ? error("bad Interval") : new(lo, hi)
    Interval(a) = Interval(a, a)
    Interval(i::Interval) = new(i.lo, i.hi)
end

# membership
function in(a::Number, x::Interval)
    x.lo <= a <= x.hi
end

# subset
function in(x::Interval, y::Interval)
    y.lo <= x.lo && x.hi <= y.hi
end

(==)(x::Interval, y::Interval) = x.lo == y.lo && x.hi == y.hi

(<)(x::Interval, y::Interval) = x.hi < y.lo
(>)(x::Interval, y::Interval) = y<x

(<)(x::RR, y::Interval) = x < y.lo
(<)(x::Interval, y::RR) = x.hi < y

(>)(x::RR, y::Interval) = y<x
(>)(x::Interval, y::RR) = y<x

zero(i::Interval) = Interval(0)
one(i::Interval) = Interval(1)

function abs(x::Interval)
    if 0 in x
        return Interval(0, max(abs(x.lo), abs(x.hi)))
    else
        return Interval(min(abs(x.lo), abs(x.hi)), max(abs(x.lo), abs(x.hi)))
    end
end

widen(lo, hi) = Interval(prevfloat(lo), nextfloat(hi))

absr(x::Interval) = Interval(RR(0), max(abs(x.lo), abs(x.hi)))

function +(x::Interval, y::Interval)
    lo::RR = x.lo+y.lo
    hi::RR = x.hi+y.hi
    widen(lo, hi)
end

(-)(x::Interval) = Interval(-x.hi, -x.lo)
(-)(x::Interval, y::Interval) = x+(-y)

function +(x::Interval, a::Number)
    lo::RR = x.lo+a
    hi::RR = x.hi+a
    widen(lo, hi)
end

(+)(a::Number, x::Interval) = x+a
(-)(x::Interval, a::Number) = x+(-a)
(-)(a::Number, x::Interval) = -x+a

function *(x::Interval, y::Interval)
    lo::RR = min(x.lo*y.lo, x.lo*y.hi, x.hi*y.lo, x.hi*y.hi)
    hi::RR = max(x.lo*y.lo, x.lo*y.hi, x.hi*y.lo, x.hi*y.hi)
    widen(lo, hi)
end

function *(x::Interval, a::Number)
    if a >= 0
        c::RR = x.lo*a
        d::RR = x.hi*a
        widen(c, d)
    else
        e::RR = x.hi*a
        f::RR = x.lo*a
        widen(e, f)
    end
end

(*)(a::Number, x::Interval) = x*a

function iRec(x::Interval)
    if 0 in x
        error("division by zero")
    else
        lo::RR = 1/x.hi
        hi::RR = 1/x.lo
        widen(lo, hi)
    end
end

(/)(x::Interval, y::Interval) = x*iRec(y)
(/)(x::Interval, a::Number) = x*(1/a)
(/)(a::Number, x::Interval) = a*iRec(x)

#handle that tightness to 1ULP is not guaranteed for Float64 ^:
@inline function sqr(x::Interval)
    if x.lo >= RR(0)
        return widen(x.lo*x.lo, x.hi*x.hi)
    elseif x.hi <= RR(0)
        return widen(x.hi*x.hi, x.lo*x.lo)
    else
        hi = nextfloat(max(x.lo*x.lo, x.hi*x.hi))
        return Interval(RR(0), hi)
    end
end

#power using the above repeated squaring
function russian_sqr(t::Interval, p::Integer)
    if p == 0
        return one(t)
    elseif p == 1
        return t
    end
    
    # ignore trailing zeros in binary expansion
    trailing = trailing_zeros(p) + 1
    p >>= trailing
    
    # repeatedly square
    while (trailing -= 1) > 0
        t = sqr(t)
    end
    
    # compute the rest
    u = t
    while p > 0
        trailing = trailing_zeros(p) + 1
        p >>= trailing
        while (trailing -= 1) >= 0
            t = sqr(t)
        end
        u *= t
    end
    u
end

#gain tight bounds by applying the above to singleton intervals
#to avoid dependency problems
function (^)(x::Interval, p::Integer)
    p == 0 && return Interval(1)
    p == 1 && return x
    p < 0 && return Interval(1)/(x^ unsigned(-p))

    if x.lo >= RR(0) || x.hi <= RR(0)
        return russian_sqr(x, p)
    end
    
    if isodd(p)
        lo = russian_sqr(Interval(x.lo), p).lo
        hi = russian_sqr(Interval(x.hi), p).hi
        return Interval(lo, hi)
    else
        hi_lo = russian_sqr(Interval(x.lo), p).hi
        hi_hi = russian_sqr(Interval(x.hi), p).hi
        return Interval(RR(0), max(hi_lo, hi_hi))
    end
end

function sqrt(x::Interval)
    @assert x.lo >= RR(0)
    lo = Base.sqrt(x.lo) #IEEE-754
    hi = Base.sqrt(x.hi) #IEEE-754
    return nonnegative(widen(lo, hi))
end

comparable(x::Interval, y::Interval) = x.hi <= y.lo || x.lo >= y.hi

(<)(x::Interval, y::Interval) = x.hi < y.lo
(>)(x::Interval, y::Interval) = y<x   

function bound_quasi_power_careful(k0, g_norm::Interval)
    @assert(g_norm.hi < Interval(1))
    @assert(g_norm.hi > Interval(0))
    # important to make the interval trivial, here
    g_norm = Interval(g_norm.hi, g_norm.hi)
    k = k0 + 1
    latest = k * g_norm ^ (k - 1)
    greatest_upper = latest.hi
    greatest_lower = latest.lo
    #println(latest)
    
    # a new upper bound strictly below a previous lower bound
    # indicates that we have passed the maximum.
    # at that point, the greatest upper bound found so far
    # and the greatest lower bound found so far
    # together bound the maximum of the graph.
    
    while latest.hi >= greatest_lower
        if latest.lo > greatest_lower
            greatest_lower = latest.lo
        end
        if latest.hi > greatest_upper
            greatest_upper = latest.hi
        end
        k = k + 1
        latest = k * g_norm ^ (k - 1)
        #println(latest)
    end
    Interval(greatest_lower, greatest_upper)
end

# here, we would make the argument that
# the rounded version of the quasi power
# is not order reversing at any point.

# this means that it should be sufficient to
# compute upper bounds only.

function bound_quasi_power(k0, g_norm::Interval)
    @assert(g_norm.hi < Interval(1))
    @assert(g_norm.hi > Interval(0))
    # important to make the interval trivial, here
    k = k0 + 1
    latest = nextfloat(k * nextfloat(g_norm.hi ^ (k - 1)))
    greatest_upper = latest
    while latest >= greatest_upper
        if latest > greatest_upper
            greatest_upper = latest
        end
        k = k + 1
        latest = nextfloat(k * nextfloat(g_norm.hi ^ (k - 1)))
    end
    Interval(0, greatest_upper)
end

function bound_quasi_power_2(c::Interval, d::Interval, M::Int)
    @assert c < RR(1) && d < RR(1)
    ceil_x_hat = Int(ceil(1/(1-c)).hi)
    biggest = d^M
    for j=2:ceil_x_hat
        k = max(0, M+1-j)
        biggest = hull(biggest, j*c^(j-1)*d^k)
    end
    biggest
end

function balanced(i::Interval)
    a = max(abs(i.lo), abs(i.hi))
    Interval(-a, a)
end
balanced(x::RR) = balanced(Interval(x))

function nonnegative(i::Interval)
    a = max(i.lo, RR(0))
    @assert i.hi>=0
    Interval(a, i.hi)
end

function hull(v::Array{Interval, 1})
    lo = min([i.lo for i in v]...)
    hi = max([i.hi for i in v]...)
    Interval(lo, hi)
end
hull(x::Interval, y::Interval) = hull([x, y])

ceil(x::Interval) = Interval(ceil(x.lo), ceil(x.hi))
    
unit_vector(t::Type, k::Integer) = [j == k ? t(1) : t(0) for j in 0:N]

(*)(A::Matrix{Interval}, c::Interval) = [x*c for x in A]
(*)(c::Interval, A::Matrix{Interval}) = A*c
