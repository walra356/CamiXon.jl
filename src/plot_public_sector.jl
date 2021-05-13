"""
    plot_matrices(data [, scale=1.0 [, select=(0,0) [; plotset="defaults"...]]])

Plot heatmap of `data`  of one or more (optionally `value-scaled`) matrices.

Keys:
* `plotset::Plot2Dset = castPlotset2D()`: object for plot parameters
* `inline::Bool = true`: switch inline/displayed input
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
function plot_matrices4(data, gain=1.0, select=(0,0);
                        plotset="defaults",
                        inline=true,
                        res=(900,600),
                        supertitle="supertitle",
                        footnote="footnote",
                        title="img",
                        textsize=10,
                        colormap=:gist_earth,
                        colorbarlabel="Intensity (counts)",
                        color_supertitle = (:black, 0.25),
                        color_footnote = (:black, 0.25))

    data = _format_matrix_array(data)                 # transforms data into a standardized array of matrices (images)

    valmin = minimum.(data)
    valmax = maximum.(data)

    dims = size.(data)                                                      # array of matrix dimensions
    plotset = plotset == "defaults" ? cast_Plotset2D.(dims) : plotset       # array of plot settings

    itrZ = select == (0,0) ? (1:length(data)) : select

    nz = length(itrZ)

    ncols = min(cld(res[1]-150,300),nz)
    nrows = inline ? cld(nz,ncols) : min(cld(nz,ncols), res[2]÷300)

    res = inline ? (res[1], max(res[2], nrows*300)) : res

    f = Figure(resolution = res)

    for i=1:nrows
        for j=1:ncols
            iz = (i-1)*ncols+j
            if iz <= nz
                set = plotset[itrZ[iz]]
                ax = Axis(f[i,j], aspect = set.aspect)
                ax.titlesize = textsize * 7/5
                ax.title = title *" $(itrZ[iz])"
                ax.xlabelsize = textsize * 6/5
                ax.ylabelsize = textsize * 6/5
                ax.xlabel = set.labels[2]
                ax.ylabel = set.labels[1]
                ax.xticklabelsize = textsize
                ax.yticklabelsize = textsize
                ax.xticks = set.ticks[2]
                ax.yticks = set.ticks[1]    # next rotate images 90 degrees to accommodate matrix convention:
                hm = heatmap!(ax, set.range[2], set.range[1], rotr90(data[itrZ[iz]])) # (1,1):(1,n) upper row
                hm.colormap = colormap
                hm.colorrange = (valmin[1], valmax[1] * set.gain * gain)
            end
        end
    end

    cb = Colorbar(f[:,ncols+1] , width = 10, limits = (0, maximum(data[1])), colormap = colormap, flipaxis = false)
    cb.ticklabelsize = textsize
    cb.labelsize = textsize * 6/5
    cb.label = colorbarlabel
    cb.flipaxis = true
    cb.height = Relative(1/2)

    strNote = " Note: output truncated at image $(ncols*nrows) (limited by specified screen resolution)"
    footnote = nz > ncols*nrows ? footnote * strNote  : footnote
    lblFoot = Label(f[nrows+1,:], footnote, textsize = textsize * 6/5, color = color_footnote)
    supertitle = Label(f[0,:], supertitle, textsize = textsize * 2 , color = color_supertitle)

    return f

end

"""
    cast_Plotset2D(dims [, gain [, aspect...]])

Set dims, gain, aspect-ratio, ticks-centering, ticks and labels of a Matrix plot.
#### Example:
```
dims = (512,30)
plot_defaults = cast_Plotset2D(dims)
  Plotset2D((512, 30), 1.0, 0.05859375, (false, false), (0.5:1.0:512.5, 0.5:1.0:30.5), (0:100:500, 0:5:30), ("y", "x"))
```
"""
function cast_Plotset2D(dims, gain=1.0, aspect=0.0, center=(false,false), steps=(0,0), labels=("y", "x"))

    typeof(aspect) <: Real || error("Aspect ratio: Real expected")
    typeof(center) == Tuple{Bool,Bool} || error("Center: Tuple{Bool,Bool} expected")
    typeof(steps) <: Tuple{Int,Int} || error("Ticks: Tuple{Int,Int} expected")
    typeof(labels) <: Tuple{String,String} || error("Labels: Tuple{String,String} expected")

    (ny,nx) = dims

    aspect = aspect == 0.0 ? nx/ny : aspect
    range = (_set_range(ny; center=center[1]),_set_range(nx; center=center[2]))
    ticks = (_set_ticks(ny, steps[1];center=center[1]), _set_ticks(nx, steps[2]; center=center[2]))

    return Plotset2D(dims, gain, aspect, center, range, ticks, labels)

end
