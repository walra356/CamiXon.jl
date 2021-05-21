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
