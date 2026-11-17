import .Base: setrounding

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
    setrounding(RR, RoundDown)
    lo::RR = x.lo+y.lo
    setrounding(RR, RoundUp)
    hi::RR = x.hi+y.hi
    setrounding(RR, RoundNearest)
    widen(lo, hi)
end

(-)(x::Interval) = Interval(-x.hi, -x.lo)
(-)(x::Interval, y::Interval) = x+(-y)

function +(x::Interval, a::Number)
    setrounding(RR, RoundDown)
    lo::RR = x.lo+a
    setrounding(RR, RoundUp)
    hi::RR = x.hi+a
    setrounding(RR, RoundNearest)
    widen(lo, hi)
end

(+)(a::Number, x::Interval) = x+a
(-)(x::Interval, a::Number) = x+(-a)
(-)(a::Number, x::Interval) = -x+a

function *(x::Interval, y::Interval)
    setrounding(RR, RoundDown)
    lo::RR = min(x.lo*y.lo, x.lo*y.hi, x.hi*y.lo, x.hi*y.hi)
    setrounding(RR, RoundUp)
    hi::RR = max(x.lo*y.lo, x.lo*y.hi, x.hi*y.lo, x.hi*y.hi)
    setrounding(RR, RoundNearest)
    widen(lo, hi)
end

function *(x::Interval, a::Number)
    if a >= 0
        setrounding(RR, RoundDown)
        c::RR = x.lo*a
        setrounding(RR, RoundUp)
        d::RR = x.hi*a
        setrounding(RR, RoundNearest)
        widen(c, d)
    else
        setrounding(RR, RoundDown)
        e::RR = x.hi*a
        setrounding(RR, RoundUp)
        f::RR = x.lo*a
        setrounding(RR, RoundNearest)
        widen(e, f)
    end
end

(*)(a::Number, x::Interval) = x*a

function iRec(x::Interval)
    if 0 in x
        error("division by zero")
    else
        setrounding(RR, RoundDown)
        lo::RR = 1/x.hi
        setrounding(RR, RoundUp)
        hi::RR = 1/x.lo
        setrounding(RR, RoundNearest)
        widen(lo, hi)
    end
end

(/)(x::Interval, y::Interval) = x*iRec(y)
(/)(x::Interval, a::Number) = x*(1/a)
(/)(a::Number, x::Interval) = a*iRec(x)

function ^(x::Interval, n::Integer)
    if n%2 == 1 || x.lo >= 0
        # odd power
        setrounding(RR, RoundDown)
        a::RR = x.lo^n
        setrounding(RR, RoundUp)
        b::RR = x.hi^n
        setrounding(RR, RoundNearest)
        return widen(a, b)
    elseif x.hi < 0
        setrounding(RR, RoundDown)
        c::RR = x.hi^n
        setrounding(RR, RoundUp)
        d::RR = x.lo^n
        setrounding(RR, RoundNearest)
        return widen(c, d)
    else
        setrounding(RR, RoundUp)
        e::RR = x.lo^n
        f::RR = x.hi^n
        setrounding(RR, RoundNearest)
        return Interval(0, nextfloat(max(e, f)))
    end
end

function sqrt(x::Interval)
    half = RR(1)/RR(2)
    @assert half*RR(2) == RR(1)
    @assert x.lo >= RR(0)
    setrounding(RR, RoundDown)
    a = x.lo^half
    setrounding(RR, RoundUp)
    b = x.hi^half
    setrounding(RR, RoundNearest)
    widen(a, b)
end

comparable(x::Interval, y::Interval) = x.hi <= y.lo || x.lo >= y.hi

# note that interval comparison needs to be handled carefully!
# if we define X<Y <=> x<y forall x,y
# and X==Y <=> x==y forall x,y (this would mean both are identical singletons)
# but we'd most likely want X<=Y <=> x<=y forall x,y
# and that is not equivalent to X<Y or X==Y.
(<)(x::Interval, y::Interval) = x.hi < y.lo
(>)(x::Interval, y::Interval) = y<x   

setrounding(t::Type{Float64}, mode) = nothing
#setrounding(t::Type{Double64}, mode) = nothing

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
    setrounding(RR, RoundUp)
    latest = nextfloat(k * nextfloat(g_norm.hi ^ (k - 1)))
    greatest_upper = latest
    while latest >= greatest_upper
        if latest > greatest_upper
            greatest_upper = latest
        end
        k = k + 1
        latest = nextfloat(k * nextfloat(g_norm.hi ^ (k - 1)))
    end
    setrounding(RR, RoundNearest)
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

# we established that the two methods produce the same results
#for k in 1:100
#    p = Int(floor(rand()*500))
#    x = 1-0.001*rand()
#    bound1 = bound_quasi_power(p, Interval(x)).hi
#    bound2 = bound_quasi_power_crude(p, Interval(x)).hi
#    if bound1 != bound2
#        println("$x: $bound1 != $bound2")
#    end
#end

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
