struct Rectangle
    re::Interval
    im::Interval
    Rectangle(re::Interval, im::Interval) = new(re, im)
    Rectangle(re, im) = new(Interval(re), Interval(im))
    Rectangle(re) = new(Interval(re), Interval(0))
end

in(x::Rectangle, z::Rectangle) = x.re in z.re && x.im in z.im
zero(z::Rectangle) = Rectangle(0, 0)
one(z::Rectangle) = Rectangle(1, 0)

(-)(z::Rectangle) = Rectangle(-z.re, -z.im)

(+)(z::Rectangle, w::Rectangle) = Rectangle(z.re+w.re, z.im+w.im)
(-)(z::Rectangle, w::Rectangle) = Rectangle(z.re-w.re, z.im-w.im)

(*)(z::Rectangle, x::Interval) = Rectangle(z.re*x, z.im*x)
(/)(z::Rectangle, x::Interval) = z*(1/x)

function reciprocal(z::Rectangle)
    d = z.re^2 + z.im^2
    Rectangle(z.re/d, -z.im/d)
end

(/)(z::Rectangle, w::Rectangle) = z*reciprocal(w)

(*)(z::Rectangle, w::Rectangle) = Rectangle(z.re*w.re-z.im*w.im, z.re*w.im+z.im*w.re)

(+)(z::Rectangle, x::Interval) = Rectangle(z.re+x, z.im)
(-)(z::Rectangle, x::Interval) = Rectangle(z.re-x, z.im)

# derived
(*)(x::Interval, z::Rectangle) = z*x
(/)(x::Interval, z::Rectangle) = Rectangle(x)/z
(+)(x::Interval, z::Rectangle) = z+x
(-)(x::Interval, z::Rectangle) = -z+x
    
#
(*)(z::Rectangle, x::RR) = z*Interval(x)
(*)(x::RR, z::Rectangle) = z*x

# we don't have sqrt of an interval yet.
abs(z::Rectangle) = sqrt(nonnegative(z.re^2 + z.im^2))

(^)(z::Rectangle, p::Integer) = russian(z, p)

Rectangle(z::Complex{BigFloat}) = Rectangle(z.re, z.im)

hull(v::Vector{Rectangle}) = Rectangle.(hull([r.re for r in v]), hull([r.im for r in v]))
