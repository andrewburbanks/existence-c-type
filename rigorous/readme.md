# Rigorous computational framework

> Implements rigorous enclosures for computations involving intervals, rectangles, and balls of functions analytic on bidisc domains.

## Important components

Preliminary:
* `Preliminaries.jl` - various operators and functions that we will overload

Nonrigorous polynomials:
* `Poly2.jl` - polynomials in two variables with truncated operations
* `Poly2D.jl` - polynomials in two variables with truncated operations, expanded with respect to a bidisc

Generic methods for both nonrigorous and rigorous versions:
* `Russian.jl` - efficient polynomial and function ball exponentiation
* `Newton.jl` - Newton method

Interval arithmetic / minimal:
* `Interval.jl`
* `Interval_Float64_crude.jl` - using Float64 and rounding to prev/next float

Interval arithmetic / additional:
* `Interval_Generic.jl` - for use with BigFloat and IEEE directed rounding modes
* `Interval_Float64_LLVM.jl` - if using Float64 with LLVM access to rounding modes
* `SetRoundingLLVM.jl`

Rectangles in the complex plane and coverings of circles, discs, bidiscs:
* `Rectangle.jl` - rectangle enclosure of complex number
* `Cover.jl` - rectangle coverings

Rigorous framework:
* `Domain2.jl` - bidisc domains in two complex variables
* `Trunc2.jl` - truncated power series with interval coefficients
* `StdBall2.jl` - standard balls of functions analytic on unit bidiscs
* `Ball2.jl` - non-standard balls of functions analytic on general bidiscs
