# =========================== Activate CairoMakie =============================

CairoMakie.activate!()

# =========================== Set Theme =======================================

theme = Theme(fontsize = 10, resolution = (900,600))
set_theme!(theme)

# ==== set_attributes(fig; title = "title", xlabel = "x", ylabel = "y") =======

function set_attributes(fig::Figure; title = "title",xlabel = "x",ylabel = "y")
    fsize = 12
    attr = (xlabelsize = 6fsize/5, ylabelsize = 6fsize/5, titlesize = 7fsize/5,
            xautolimitmargin = (.09,.09), yautolimitmargin = (.09,.09), 
            titlefont = "TeX Gyre Heros Makie Normal, ")
    return (attr..., title = title, xlabel = xlabel, ylabel = ylabel, )
end

println("CairoMakie.activate!()")
println("theme = Theme(fontsize = 10, resolution = (900,600))")
println("set_theme!(theme)")
println("set_attributes: 
fsize = 12, 
attr = (xlabelsize = 6fsize/5, ylabelsize = 6fsize/5, titlesize = 7fsize/5,
            xautolimitmargin = (.1,.1), yautolimitmargin = (.1,.1), )")

