function _auto_ticks(n::Real, center::Bool)

    if 1 <= n <= 10
         ticks = 1
    elseif 10 < n <=14
         ticks = 2
    elseif 14 < n <=20
         ticks = 5
    else
        n = center ? n÷2 : n
        p = exp10(log10(n)-floor(Int,log10(n)))
        s = 10^(length(string(n))-1)
        if p > 8.4
            ticks = center ? 4s : 2s
        elseif p > 6.3
            ticks = center ? 3s : s
        elseif p > 4.2
            ticks = center ? 2s : s
        elseif p > 2.1
            ticks = center ? s : s÷2
        else
            ticks = center ? 5(s÷10) : 2(s÷10)
        end
    end

    return ticks

end

_irow(i::Int,ncols::Int) = (i-1)÷ncols+1

_icol(i::Int,ncols::Int) = Base.mod1(i, ncols)

_range(n::Real, center=false) = center ? (Base.isodd(n) ? (-n/2:1:n/2) : (0.5-n/2:1:0.5+n/2)) : (0.5:1:0.5+n)

function _ticks(n::Real, step=0, center=false)

    step = Base.iszero(step) ? _auto_ticks(n, center) : step

    max = (n÷2÷step)step

    return center ? (-max:step:max) : (0:step:n+1)

end
