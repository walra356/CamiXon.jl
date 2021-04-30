"""
    plot_matrix(data [, scale=1.0 [; aspect=0...]])

Keys:
* `aspect::Real=0`: aspect ratio (default: not scaled)
* `inline::Bool=true`: switch between inline and displayed input (default: inline)
* `res::Tuple`: screen resolution (default: '(900,600)')
* `note::String`: footnote text
* `textsize::Int`: textsize (default: 10pt)
* `center::Tuple{Bool,2}`: centering of scale on (x axis, y axis)
#### Example:
```
```
"""
function plot_matrix(data, scale=1; aspect=0, inline=true, res=(900,600), note=" ", textsize=10, center=(false,false))

    typeof(data) <: Matrix || error("'$data': is not a 2D image")

    data = rotr90(data)

    GLMakie.activate!()

    AbstractPlotting.inline!(inline)

    valmin = minimum(data)
    valmax = maximum(data) * scale

    nx,ny = size(data)

    aspect = aspect == 0 ? nx/ny : 1

    (xticks,yticks) = _set_ticks(data; center=(center[1],center[2]))

    itrX = center[1] ? (1-nx÷2:nx÷2) : (0:nx)
    itrY = center[2] ? (1-ny÷2:ny÷2) : (0:ny)

    res = inline ? (res[1],max(res[2],300)) : res

    (scene, layout) = layoutscene(resolution = res)

    ax = layout[1,1] = Axis(scene, aspect=aspect)
    ax.titlesize = textsize * 7/5
    ax.xlabelsize = textsize * 6/5
    ax.ylabelsize = textsize * 6/5
    ax.xlabel = "x"
    ax.ylabel = "y"
    ax.xticklabelsize = textsize
    ax.yticklabelsize = textsize
    ax.xticks = xticks
    ax.yticks = yticks

    hm = heatmap!(ax, itrX, itrY, data)
    hm.colormap = :gist_earth
    hm.colorrange = (valmin,valmax)

    cb = layout[1,2] = Colorbar(scene, width = 20, limits=(valmin,valmax), colormap = hm.colormap, flipaxis = false)
    cb.ticklabelsize = textsize
    cb.labelsize = textsize * 6/5
    cb.label = "value"
    cb.flipaxis = true

    lblFoot = layout[2,:] = Label(scene, note, textsize = textsize * 6/5)
    supertitle = layout[0,:] = Label(scene, "supertitle", textsize = textsize * 2 , color = (:black, 0.25))

    return inline ? scene : display(scene)

end
