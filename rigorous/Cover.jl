function cover_unit_square(n::Int)
    imag_unit = Rectangle(Interval(0), Interval(1))
    vals = LinRange(-1, 1, n+1)
    rects = Rectangle[]
    for j=1:length(vals)-1
        for k=1:length(vals)-1
            X = Interval(vals[j], vals[j+1])
            Y = Interval(vals[k], vals[k+1])
            R = X + Y*imag_unit
            push!(rects, R)
        end
    end
    rects
end

function cover_unit_disc(n::Int)
    rects = Rectangle[]
    for R in cover_unit_square(n)
        if ~((R.re)^2 + (R.im)^2 > RR(1))
            push!(rects, R)
        end
    end
    rects
end

function cover_unit_circle(n::Int)
    rects = Rectangle[]
    for R in cover_unit_square(n)
        if RR(1) in (R.re)^2 + (R.im)^2
            push!(rects, R)
        end
    end
    rects
end
