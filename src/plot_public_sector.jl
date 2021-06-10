"""
    plot_matrices(data [, ncols=3 [, select=(0,0) [; kwargs...]]])

Plot heatmap of `data`  of one or more (optionally `value-scaled`) matrices.

Keyword arguments:
* `gain::Real = 1.0`: for scaling the value of the matrix elements
* `steps::Tuple{Real,Real} = (0,0)`: x and y ticks stepsizes
* `center::Tuple{Bool,Bool} = (false,false)`: x and y centering of ticks ranges
* `res::Tuple{Int,Int} = (900,600)`: screen resolution
* `supertitle::String = "supertitle"`
* `color_supertitle::Tuple{Symbol,Real} = (:black, 0.25)`
* `footnote::String = "footnote"`
* `color_footnote::Tuple{Symbol,Real} = (:black, 0.25)`
* `title::String = "img"`
* `textsize::Int = 10`
* `colormap::Symbol = :gist_earth`
* `colorbarlabel::String = "Intensity (units)"`: label for scale of colorbar
"""
function plot_matrices(σ, ncols=3, select=(0,0);
                        gain=1.0,
                        steps=(1,1),
                        center=(false,false),
                        res=(900,600),
                        supertitle="supertitle",
                        footnote="footnote",
                        title="img",
                        textsize=10,
                        colormap=:gist_earth,
                        colorbarlabel="Intensity (counts)",
                        color_supertitle = (:black, 0.4),
                        color_footnote = (:black, 0.4))


    σ = rotr90.(σ)                   # rotate matrix to obtain a template image of the matrix

    valmin = minimum.(σ)
    valmax = maximum.(σ)

    itrZ = select == (0,0) ? (1:length(σ)) : select
    dims = size.(σ)                   # array of matrix dimensions

    nx = [dims[i][1] for i ∈ itrZ]
    ny = [dims[i][2] for i ∈ itrZ]

    nz = length(itrZ)

    nrows = _irows(nz,ncols)

    fig = Figure(resolution=res)

    ax = [Axis(fig[_irows(i,ncols),_icols(i,ncols)]) for i=1:nz]
    hm = [heatmap!(ax[i], _range(nx[i],center[1]), _range(ny[i],center[2]),σ[itrZ[i]]) for i ∈ eachindex(ax)]

    for i ∈ eachindex(ax)
        ax[i].aspect = nx[i]/ny[i]
        ax[i].titlesize = textsize * 7/5
        ax[i].title = title *" $(itrZ[i])"
        ax[i].xlabelsize = textsize * 6/5
        ax[i].ylabelsize = textsize * 6/5
        ax[i].xlabel = "x"
        ax[i].ylabel = "y"
        ax[i].xticklabelsize = textsize
        ax[i].yticklabelsize = textsize
        ax[i].xticks = _ticks(nx[i], steps[1], center[1])
        ax[i].yticks = _ticks(ny[i], steps[2], center[2])
        hm[i].colormap = colormap
        hm[i].colorrange = (valmin[1], valmax[1] * gain)
    end

    cb = Colorbar(fig[:,ncols+1] , width = 10, limits = (0, valmax[1]), colormap = colormap, flipaxis = false)
    cb.ticklabelsize = textsize
    cb.labelsize = textsize * 6/5
    cb.label = colorbarlabel
    cb.flipaxis = true
    cb.height = Relative(1/2)

    strNote = " Note: output truncated at image $(ncols*nrows) (limited by specified screen resolution)"
    footnote = nz > ncols*nrows ? footnote * strNote  : footnote
    lblFoot = Label(fig[nrows+1,:], footnote, textsize = textsize * 6/5, color = color_footnote)
    supertitle = Label(fig[0,:], supertitle, textsize = textsize * 2 , color = color_supertitle)

    return fig

end


# ==================================== step125(x) ============================================================

"""
    step125(x)

Step used for deviding the number x in steps according to 1-2-5 scheme
#### Examples:
```
step125.([5,10,21.3,50,100.1])
5-element Vector{Int64}:
  1
  2
  5
 10
 20
```
"""
function step125(x::Real)

    m = log10_mantissa(x)
    p = log10_characteristic_power(x)
    v = 10^m
    d = v > 7.9 ? 2.0 : v > 3.9 ? 1.0 : v > 1.49 ? 0.5 : 0.2

    return max(1,round(Int, d *= 10^p))

end


# ==================================== select125(x) ===============================================================

"""
    select125(x)

Select elements of the collection x by index according to 1-2-5 scheme
#### Examples:
```@docs
x1 = [1,2,4,6,8,10,13,16,18,20,40,60,80,100]
x2 = string.(x1)
x3 = Base.OneTo(100)
x4 = (1:100)
println(select125(x1)); println(select125(x2)); println(select125(x3)); println(select125(x4))
 [2, 6, 10, 16, 20, 60, 100]
 ["2", "6", "10", "16", "20", "60", "100"]
 [20, 40, 60, 80, 100]
 [20, 40, 60, 80, 100]
```
"""
select125(x) = (n = length(x); return [x[i] for i=step125(n):step125(n):n])


# ==================================== ticks125(x) ================================================================

"""
    ticks125(x)

Tickvalues for x according to 1-2-5 scheme
#### Examples:
```
ticks125.([5,10,21.3,50,100.1])
5-element Vector{StepRange{Int64, Int64}}:
 -5:1:5
 -10:2:10
 -20:5:20
 -50:10:50
 -100:20:100
```
"""
function ticks125(x)

    Δx = step125(x)
    xm = (round(Int, x) ÷ Δx) * Δx
    return (-xm:Δx:xm)

end


# ==================================== centerticks(pos) ============================================================

"""
    centerticks(pos)

Tick positions centered with respect to the pixels of a heatmap (1-2-5 scheme)
Note that ticks outside the plot limits are not rendered.
#### Examples:
```
x1 = Base.OneTo(100)
x2 = UnitRange(-200:100)
x3 = -200..100
x4 = [1,4,2,5] # specification of segment lenghts
println(centerticks(x1)); println(centerticks(x2)); println(centerticks(x3)); println(centerticks(x4))
 -100:20:100
 -200:50:200
 -200:50:200
 Float32[0.5, 3.0, 6.0, 9.5]
```
"""
centerticks(pos::Base.OneTo{Int64}) = ticks125(maximum(abs.(collect(pos))))
centerticks(pos::UnitRange{Int}) = ticks125(maximum(abs.(collect(pos))))
centerticks(pos::IntervalSets.ClosedInterval{<:Real}) = ticks125(maximum(abs.(collect(pos.left:pos.right))))
centerticks(seg::Vector{<:Real}) =  select125([Base.sum(seg[1:i]) for i ∈ Base.eachindex(seg)] .- seg .* 0.5f0)


# ==================================== edgeticks(pos) ============================================================

"""
    edgeticks(pos)

Tick positions corresponding to the pixel edges of a heatmap (1-2-5 scheme)
Note that ticks outside the plot limits are not rendered.
#### Examples:
```
x1 = Base.OneTo(100)
x2 = UnitRange(-200:100)
x3 = -200..100
x4 = [1,4,2,5] # specification of segment lenghts
println(edgeticks(x1)); println(edgeticks(x2)); println(edgeticks(x3)); println(edgeticks(x4))
 -100:20:100
 -200:50:200
 -200:50:200
 [0, 1, 5, 7, 12]
```
"""
edgeticks(pos::Base.OneTo{Int64}) = centerticks(pos)
edgeticks(pos::UnitRange{Int}) = centerticks(pos)
edgeticks(pos::IntervalSets.ClosedInterval{<:Real}) = centerticks(pos)
edgeticks(seg::Vector{<:Real}) = Base.append!([0],[Base.sum(seg[1:i]) for i ∈ Base.eachindex(seg)])


# ==================================== centers(pos) ============================================================

"""
    centers(pos)

Positions centered with respect to the pixels of a heatmap (1-2-5 scheme).
#### Examples:
```
x1 = Base.OneTo(100)
x2 = UnitRange(-200:100)
x3 = -200..100
x4 = [1,4,2,5] # specification of segment lenghts
println(edgeticks(x1)); println(edgeticks(x2)); println(edgeticks(x3)); println(edgeticks(x4))
 Base.OneTo(100)
 -200:100
 -200..100
 [0, 1, 5, 7, 12]
```
"""
centers(pos::Base.OneTo{Int64}) = pos
centers(pos::UnitRange{Int}) = pos
centers(pos::IntervalSets.ClosedInterval{<:Real}) = pos
centers(seg::Vector{<:Real}) = (s=Base.append!([0],seg); return [Base.sum(s[1:i]) for i ∈ Base.eachindex(s)])

"""
    edges(pos)

Positions corresponding to the pixel edges of a heatmap (1-2-5 scheme).
#### Examples:
```
x1 = Base.OneTo(100)
x2 = UnitRange(-200:100)
x3 = -200..100
x4 = [1,4,2,5] # specification of segment lenghts
println(edges(x1)); println(edges(x2)); println(edges(x3)); println(edges(x4))
 0.5..99.5
 0.5..99.5
 -200.5..99.5
 [0, 1, 5, 7, 12]
```
"""
edges(pos::Base.OneTo{Int64}) = (0.5)..(pos[end] - 0.5)
edges(pos::UnitRange{Int}) = (0.5)..(pos[end] - 0.5)
edges(pos::IntervalSets.ClosedInterval{<:Real}) = (pos.left-0.5)..(pos.right-0.5)
edges(seg::Vector{<:Real}) = (s=Base.append!([0],seg); return [Base.sum(s[1:i]) for i ∈ Base.eachindex(s)])
