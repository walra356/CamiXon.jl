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

function _format_matrix_array(data)

    datatype = typeof(data)

    ismatrix_array = datatype <: Array{Matrix{T}, 1} where T<:Real ? true : false
    isarray1D = datatype <: Array{T,1} where T<:Real ? true : false
    isarray2D = datatype <: Array{T,2} where T<:Real ? true : false
    isarray3D = datatype <: Array{T,3} where T<:Real ? true : false
    isarray = datatype <: (Array{T,N} where {T,N}) ? true : false

    isarray || error("PlotError: datatype not a 1D, 2D or 3D array")
    data = isarray1D ? (ismatrix_array ? data : error("PlotError: datatype not a matrix array")) : data
    data = isarray2D ? [data] : data
    data = isarray3D ? [data[:,:,i] for i=1:size(data)[3]] : data

    return data

end

function _set_range(n::Real, center=false)

    return center ? (isodd(n) ? (-n/2:1:n/2) : (0.5-n/2:1:0.5+n/2)) : (0.5:1:0.5+n)  # return range

end

function _set_ticks(n::Real, center=false, step=0)

    step = step == 0 ? _auto_ticks(n, center) : step

    max = (n÷2÷step)step

    return center ? (-max:xticks:max) : (0:step:n)  # return ticks

end
