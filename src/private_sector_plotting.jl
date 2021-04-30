function _auto_ticks(n::Int; center=false)

    if 1 <= n <= 10
         ticks = center ? (-n÷2:1:n÷2) : (0:1:n)
    elseif 10 < n <=14
         ticks = center ? (-n÷2:2:n÷2) : (0:2:n)
    elseif 14 < n <=20
         ticks = center ? (-n÷2:5:n÷2) : (0:5:n)
    else
        n = center ? n÷2 : n
        p = exp10(log10(n)-floor(Int,log10(n)))
        s = 10^(length(string(n))-1)
        if p > 8.4
            ticks = center ? (-8:4:8).*s : (0:2s:n)
        elseif p > 6.3
            ticks = center ? (-6:3:6).*s : (0:s:n)
        elseif p > 4.2
            ticks = center ? (-4:2:4).*s : (0:s:n)
        elseif p > 2.1
            ticks = center ? (-2:1:2).*s : (0:(s÷2):n)
        else
            ticks = center ? (-10:5:10).*(s÷10) : (0:2(s÷10):n)
        end
    end

    return ticks

end

function _set_ticks(m::Matrix; center=(false,false), xticks=(1:1:1), yticks=(1:1:1))

    (ny,nx) = size(m)

    xticks = xticks == (1:1:1) ? _auto_ticks(nx; center=center[1]) : xticks
    yticks = yticks == (1:1:1) ? _auto_ticks(ny; center=center[2]) : yticks

    return (yticks,xticks)

end
