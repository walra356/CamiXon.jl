function _test_saturation(data; bitchip::Int=14, bitpix::Int=32, bzero::Int=0, numkin::Int=1)

    bitchip_sat::Int = (2^bitchip - 1) * numkin
    bitpix_sat::Int = 2^bitpix - bzero

    data = convert(Array{Int,3},data)

    o1::Array{Int,1} = []
    o2::Array{Int,1} = []

    [append!(o1, CamiXon.find_all(vec(data[:,:,i]), bitchip_sat; count=true)[1] for i ∈ axes(data,3))]
    [append!(o2, CamiXon.find_all(vec(data[:,:,i]), bitpix_sat; count=true)[1] for i ∈ axes(data,3))]

    sum(o1) > 0 ? println("Warning: physical saturation detected ($(sum(o1)) times)") : false
    sum(o2) > 0 ? println("Warning: $bitpix bit saturation detected ($(sum(o2)) times)") : false

end


function _fits_scene(filnam::String, img=(0,0); inline=true, res=(900,600), footnote::String=" ", textsize=10, xcenter=false, ycenter=false)

    GLMakie.activate!()

    AbstractPlotting.inline!(inline)

    isfile(filnam) ? true : error(filnam * ": not found in current directory")

    f = fits_read(filnam)
    data = f[1].dataobject.data
    dict = f[1].header.dict

    bitchip = 14
    colormap = :gist_earth

    nx = get(dict,"NAXIS1",0)
    ny = get(dict,"NAXIS2",0)
    nz = get(dict,"NAXIS3",0)
    bp = get(dict,"BITPIX",0)
    bz = get(dict,"BZERO",0)
    nk = get(dict,"NUMKIN",0)

    val_sat = (2^bitchip - 1) * nk

    _test_saturation(data; bitchip=bitchip, bitpix=bp, bzero=bz, numkin=nk)

    (yticks,xticks) = _set_ticks(data[:,:,1]; xcenter=xcenter, ycenter=ycenter)

    itrX = xcenter ? (1-nx÷2:nx÷2) : (0:nx)
    itrY = ycenter ? (1-ny÷2:ny÷2) : (0:ny)

    itrZ = img == (0,0) ? (1:nz) : img
    nz = length(itrZ)

    ncols = min(res[1]÷266,nz)
    nrows = inline ? cld(nz,ncols) : min(cld(nz,ncols),res[2]÷300)

    res = inline ? (res[1],max(res[2],nrows*300)) : res

    (scene, layout) = layoutscene(resolution = res)

    for i=1:nrows
        for j=1:ncols
            iz = (i-1)*ncols+j
            if iz <= nz
                ax = layout[i,j] = Axis(scene, aspect = 1)
                ax.titlesize = textsize * 7/5
                ax.title = "img $(itrZ[iz])"
                ax.xlabelsize = textsize * 6/5
                ax.ylabelsize = textsize * 6/5
                ax.xlabel = "x (pixel)"
                ax.ylabel = "y (pixel)"
                ax.xticklabelsize = textsize
                ax.yticklabelsize = textsize
                ax.xticks = xticks
                ax.yticks = yticks
                hm = heatmap!(ax, itrX, itrY, data[:,:,itrZ[iz]])
                hm.colormap = colormap
            end
        end
        cb = layout[i,ncols+1] = Colorbar(scene, width = 10, limits = (0, val_sat), colormap = colormap, flipaxis = false)
        cb.ticklabelsize = textsize
        cb.labelsize = textsize * 6/5
        cb.label = "intensity (counts)"
        cb.flipaxis = true
    end

    strNote = " Note: output truncated at image $(ncols*nrows) (limited by specified screen resolution)"
    footnote = nz > ncols*nrows ? footnote * strNote  : footnote
    lblFoot = layout[nrows+1,:] = Label(scene, footnote, textsize = textsize * 6/5)
    supertitle = layout[0,:] = Label(scene, filnam, textsize = textsize * 2 , color = (:black, 0.25))

    return scene

end

function fits_display222(filnam::String, img=(0,0); res=(900,600), footnote::String=" ", textsize=10, xcenter=false, ycenter=false)

    GLMakie.activate!()

    sc = _fits_scene(filnam, img; inline=false, res=res, footnote=footnote, textsize=textsize, xcenter=xcenter, ycenter=ycenter)

    return display(sc)

end

"""
    plot_matrix(data [, scale=1.0 [; aspect=0...]])

Heatmap of the (optionally `scaled`) matrix `data` using the colormap `:gist_earth`.

Keys:
* `aspect::Real=0`: aspect ratio (default: not scaled)
* `inline::Bool=true`: switch between inline and displayed input (default: inline)
* `res::Tuple`: screen resolution (default: '(900,600)')
* `note::String`: footnote text
* `textsize::Int`: textsize (default: 10pt)
* `center::Tuple{Bool,2}`: centering of scale on (x axis, y axis)
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
