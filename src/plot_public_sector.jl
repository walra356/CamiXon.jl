const PLOT_DEFAULTS = cast_Plot2Dset()

"""
    plot_matrices(data [, scale=1.0 [, select=(0,0) [; plotset= [, inline= [, supertitle= [, footnote= [, res=, [testsize=, [colormap=]]]]]]]]])

Heatmap of `data`  of one or more (optionally `value-scaled`) matrices.

Keys:
* `plotset::Plot2Dset=PLOT_DEFAULTS`: aspect=0, center=(true, true), ticks=(1:1:1, 1:1:1), title="title", labels=("x", "y")
* `inline::Bool=true`: switch between inline and displayed input (default: `inline`)
* `supertitle::String`: supertitle (default = `"supertitle"`)
* `footnote::String`: footnote (default = `"footnote"`)
* `res::Tuple`: screen resolution (default: `(900,600)`)
* `textsize::Int`: textsize (default: `10`)
* `colormap::Symbol`: colormap (default: `:gist_earth`)
"""
function plot_matrices(data, scale=1, select=(0,0); plotset=PLOT_DEFAULTS, supertitle="supertitle", footnote="footnote", settings=PLOT_DEFAULTS, inline=true, res=(900,600))

    CairoMakie.activate!()

    AbstractPlotting.inline!(inline)

    data = _format_matrix_array(data)                 # transforms data into a standardized array of matrices
    data = rotr90.(data)

    set = plotset == PLOT_DEFAULTS ? PLOT_DEFAULTS : set

    (nx,ny) = size(data[1])

    aspect = set.aspect == 0 ? nx/ny : set.aspect

    itrX = set.center[1]   ? (1-nx÷2:nx÷2)    : (0:nx)
    itrY = set.center[2]   ? (1-ny÷2:ny÷2)    : (0:ny)
    itrZ = select == (0,0) ? (1:length(data)) : select

    nz = length(itrZ)

    ncols = min(res[1]÷266,nz)
    nrows = inline ? cld(nz,ncols) : min(cld(nz,ncols),res[2]÷300)

    res = inline ? (res[1], max(res[2], nrows*300)) : res

    (scene, layout) = layoutscene(resolution = res)
    (xticks,yticks) = _set_ticks(data[1]; set.center, set.ticks)
    (valmin,valmax) = (minimum(data[1]), maximum(data[1]) * scale)

    for i=1:nrows
        for j=1:ncols
            iz = (i-1)*ncols+j
            if iz <= nz
                ax = layout[i,j] = Axis(scene, aspect = aspect)
                ax.titlesize = textsize * 7/5
                ax.title = "img $(itrZ[iz])"
                ax.xlabelsize = textsize * 6/5
                ax.ylabelsize = textsize * 6/5
                ax.xlabel = set.labels[1]
                ax.ylabel = set.labels[2]
                ax.xticklabelsize = textsize
                ax.yticklabelsize = textsize
                ax.xticks = xticks
                ax.yticks = yticks
                hm = heatmap!(ax, itrX, itrY, data[itrZ[iz]])
                hm.colormap = colormap
                hm.colorrange = (valmin,valmax)
            end
        end
        cb = layout[i,ncols+1] = Colorbar(scene, width = 10, limits = (0, valmax), colormap = colormap, flipaxis = false)
        cb.ticklabelsize = textsize
        cb.labelsize = textsize * 6/5
        cb.label = "magnitude"
        cb.flipaxis = true
    end

    strNote = " Note: output truncated at image $(ncols*nrows) (limited by specified screen resolution)"
    footnote = nz > ncols*nrows ? footnote * strNote  : footnote
    lblFoot = layout[nrows+1,:] = Label(scene, footnote, textsize = set.textsize * 6/5)
    supertitle = layout[0,:] = Label(scene, "supertitle", textsize = set.textsize * 2 , color = (:black, 0.25))

    return inline ? scene : display(scene)

end

function cast_Plot2Dset(aspect=0, center=(true, true), ticks=(1:1:1, 1:1:1), title="title", labels=("x", "y"))

    typeof(aspect) <: Real || error("Aspect ratio: Real expected")
    typeof(center) == Tuple{Bool,Bool} || error("Center: Tuple{Bool,Bool} expected")
    typeof(ticks) <: Tuple{AbstractRange,AbstractRange} || error("Ticks: Tuple{AbstractRange,AbstractRange} expected")
    typeof(title) == String || error("Title: String expected")
    typeof(labels) <: Tuple{String,String} || error("Labels: Tuple{String,String} expected")

    return Plot2Dset(aspect, center, ticks, title, labels)

end
