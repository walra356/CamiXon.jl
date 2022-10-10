# =============================================================================
function plot_tabulated_function(Nstart::Int, Nstop::Int, f::Vector{T}, grid::Grid{T}; title="", theme = Theme(fontsize = 10, resolution = (900,350)))  where T <: Real

    set_theme!(theme)

        r = grid.r
        N = grid.N
        n = [i for i=1:N]

    1 ≤ Nstart ≤ N || error("Error: Nstart outside index range 1:$(N)")
    1 ≤ Nstop  ≤ N || error("Error: Nstop outside index range 1:$(N)")

    ylabel1 = ("f(r)") 
    ylabel2 = ("f(n)")

    itr0 = Nstart:Nstop

    fig = Figure()

    attr1 = set_attributes(fig; title = ylabel1 * " on Grid[$(itr0)]", xlabel = "r (a.u.)", ylabel = ylabel1)
    attr2 = set_attributes(fig; title = ylabel2 * " on Grid[$(itr0)]", xlabel = "n", ylabel = ylabel2)

    ax0 = Label(fig; text = title, textsize = 24,  color=:gray)
    ax1 = Axis(fig; attr1...)
    ax2 = Axis(fig; attr2...)

    lines!(ax1, r[itr0], f[itr0], markersize = 1, color=:gray90)
    lines!(ax2, n[itr0], f[itr0], markersize = 1, color=:gray90)

    lines!(ax1, r[itr0], f[itr0], markersize = 1, color=:gray90)
    scatter!(ax1, r[itr0], f[itr0], markersize = 2)

    lines!(ax2, n[itr0], f[itr0], markersize = 1, color=:gray90)
    scatter!(ax2, n[itr0], f[itr0], markersize = 2)

    fig[1,1] = ax1
    fig[1,2] = ax2
    fig[0,:] = ax0

    return fig

end
function plot_tabulated_function(itr::UnitRange{Int}, f::Vector{T}, grid::Grid{T}; title=nothing, theme = Theme(fontsize = 10, resolution = (900,350)))  where T <: Real

    Nstart = itr.start
    Nstop  = itr.stop

    return plot_tabulated_function(Nstart, Nstop, f, grid; title, theme)

end

# ===================== plot_gridfunction(itr, grid; title=nothing) ===========

function plot_gridfunction(itr::UnitRange{Int}, grid::Grid{T}; theme = Theme(fontsize = 10, resolution = (900,600)), title=nothing) where T <: Number

    set_theme!(theme)    

    N = grid.N
    1 ≤ itr.start ≤ N || error("Error: Nstart outside index range 1:$(N)")
    1 ≤ itr.stop  ≤ N || error("Error: Nstop outside index range 1:$(N)")

    # strID = length(title) > 0 ? title : grid.name * " grid"
    strID = isnothing(title) ? grid.name * " grid" : title

    r = grid.r
    r′= grid.r′

    fig = Figure()                 # start plot

    attr = set_attributes(fig)

    ax0 = Label(fig; text = strID, textsize = 24,  color=:gray)
    ax1 = Axis(fig; attr..., title = "radial distance (log scale)", xlabel = "n", ylabel = "r (a.u.)", yscale=log10)
    ax2 = Axis(fig; attr..., title = "radial derivative (log scale)", xlabel = "n", ylabel = "dr/dn (a.u.)", yscale=log10)
    ax3 = Axis(fig; attr..., title = "radial distance (lin scale)", xlabel = "n", ylabel = "r (a.u.)")
    ax4 = Axis(fig; attr..., title = "radial derivative (lin scale)", xlabel = "n", ylabel = "dr/dn (a.u.)")

    scatter!(ax1, itr, r[itr] .+ 0.000001, markersize = 2)      # avoid crash on log scale
    scatter!(ax2, itr, r′[itr] .+ 0.000001, markersize = 2)
    scatter!(ax3, itr, r[itr],  markersize = 2)
    scatter!(ax4, itr, r′[itr],  markersize = 2)

    fig[1,1] = ax1
    fig[1,2] = ax2
    fig[2,1] = ax3
    fig[2,2] = ax4
    fig[0,:] = ax0

    return fig

end
function plot_gridfunction(Nstart::Int, Nstop::Int, grid::Grid{T}; title=nothing) where T <: Number

    itr = Nstart:Nstop

    return plot_gridfunction(itr, grid; title)

end

# =============================================================================

function find_energy_interval(Emin, Emax, def::Def{T}) where T <: Number
# =============================================================================
#  find indices of energy interval on potential curve
# =============================================================================
    Nb = def.pos.Nb
    v = def.pot
    s = def.scr

    pot = v .+ s

    Emin = convert(T, Emin)
    Emax = convert(T, Emax)
    Eref = convert(T, -5)       # reference for "clipping" of amplitude

    n = Nb
    while Emin < pot[n]
        n > 1 ? n -= 1 : break
    end

    Nmax = n

    while Eref < pot[n] < Emax
        n > 1 ? n -= 1 : break
    end

    Nmin = n ≠ Nmax ? n : Nmax÷2

    return Nmin, Nmax

end
# .............................................................................
function plot_potentials(E::T, grid::Grid{T}, def::Def{T}; theme = Theme(fontsize = 10, resolution = (900,600))) where T <: Real

    set_theme!(theme)    

         N = grid.N
         r = grid.r
    symbol = def.atom.element.symbol
      name = def.orbit.name
         ℓ = def.orbit.ℓ
         v = def.pot
         s = def.scr

    pot = v .+ s

    Nlctp = get_Nlctp(E, def)
     Nmin = get_Nmin(def)
    Nuctp = get_Nuctp(E, def)

    Nstart, Nstop = find_energy_interval(5E, 0.0, def)

    itr = Nstop:N

    # println("Energy interval plot: ", itr)

    fig = Figure()                   # start plot

    attr = set_attributes(fig)

    ax0 = Label(fig; text = "$(symbol) $(name)", textsize = 24,  color=:gray)
    ax1 = Axis(fig; attr..., title = "V effective ", xlabel = "n", ylabel = "Veff (Hartree)")
    ax2 = Axis(fig; attr..., title = "V screening", xlabel = "n", ylabel = "Vscr (Hartree)")
    ax3 = Axis(fig; attr..., title = "V total", xlabel = "n", ylabel = "Vtot (Hartree)")
    ax4 = Axis(fig; attr..., title = "V total", xlabel = "r (a.u.)", ylabel = "Vtot (Hartree)")
    ax5 = Axis(fig; attr..., title = "V total", xlabel = "n", ylabel = "Vtot (Hartree)")

    lines!(ax1, itr, v[itr], markersize = 1)
    scatter!(ax1, itr, v[itr], markersize = 3)
    lines!(ax1, [Nstop,N], [0,0], markersize = 1, color=:Grey)
    lines!(ax1, [Nstop,N], [E,E], markersize = 1, color=:Green)
    E < pot[end] ? scatter!(ax1, (Nuctp, E), markersize = 5, color=:Black) : false

    lines!(ax2, itr, s[itr], markersize = 1)
    scatter!(ax2, itr, s[itr], markersize = 3)
    lines!(ax2, [Nstop,N], [0,0], markersize = 1, color=:Grey)
    lines!(ax2, [Nstop,N], [E,E], markersize = 1, color=:Green)

    lines!(ax3, itr, pot[itr], markersize = 1)
    scatter!(ax3, itr, pot[itr], markersize = 3)
    lines!(ax3, [Nstop,N], [0,0], markersize = 1, color=:Grey)
    lines!(ax3, [Nstop,N], [E,E], markersize = 1, color=:Green)
    E < pot[end] ? scatter!(ax3, (Nuctp, E), markersize = 5, color=:Black) : false

    ℓ > 0 ? Emin = minimum(pot) : Emin = -0.1

    Nstart, Nstop = find_energy_interval(Emin,-0.1Emin, def)
    itr = Nstart:N

    # println("Energy interval plot: ", itr)

    lines!(ax4, r[itr], pot[itr], markersize = 1)
    scatter!(ax4, r[itr], pot[itr], markersize = 3)
    lines!(ax4, [r[1],r[N]], [0,0], markersize = 1, color=:Grey)
    lines!(ax4, [r[1],r[N]], [E,E], markersize = 1, color=:Green)
    E < pot[end] ? scatter!(ax4, (r[Nuctp], E), markersize = 5, color=:Black) : false

    lines!(ax5, itr, pot[itr], markersize = 1)
    scatter!(ax5, itr, pot[itr], markersize = 3)
    lines!(ax5, [1,N], [0,0], markersize = 1, color=:Grey)
    lines!(ax5, [1,N], [E,E], markersize = 1, color=:Green)
    E < pot[end] ? scatter!(ax5, (Nuctp, E), markersize = 5, color=:Black) : false

    fig[1,1] = ax1
    fig[1,2:3] = ax2
    fig[1,4] = ax3
    fig[2,1:2] = ax4
    fig[2,3:4] = ax5
    fig[0,:] = ax0

    str = string(r[Nuctp])
    lst = min(13, length(str))

    println("Nlctp = $(Nlctp), Nmin = $(Nmin), Nuctp = $(Nuctp) (Ructp = " * str[1:lst] * " a.u.)")

    return fig

end

# =============================================================================

function _zone(n::Int, def::Def)

        k = def.k
       Na = def.pos.Na
    Nuctp = def.pos.Nuctp
       Nb = def.pos.Nb
        N = def.pos.N

        0 < n ≤ k+1    && return 1
        k < n ≤ Na     && return 2
       Na < n ≤ Nuctp  && return 3
    Nuctp < n ≤ Nb     && return 4
       Nb < n ≤ N-k-1  && return 5
    N-k-1 < n ≤ N      && return 6

    error("Error: index outside of range 1:$(N)")

end

function _zone_color(zone::Int)

    zone == 1 && return :Black
    zone == 2 && return :Black
    zone == 3 && return :Royalblue
    zone == 4 && return :Red
    zone == 5 && return :Indianred
    zone == 6 && return :Black

    error("Error: color zone index outside of range 1:6")

end

function _zone_markersize(zone::Int)

    zone == 1 && return 2.0
    zone == 2 && return 2.0
    zone == 3 && return 1.5
    zone == 4 && return 1.5
    zone == 5 && return 1.5
    zone == 6 && return 2.0

    error("Error: color zone index outside of range 1:6")

end

function _scatterX!(ax, r, X, itr, def)

        k = def.k
       Na = def.pos.Na
    Nuctp = def.pos.Nuctp
       Nb = def.pos.Nb
        N = def.pos.N

    N1 = itr.start
    N2 = itr.stop

    A = _zone(N1, def)
    B = _zone(N2, def)

    if B == A
        itrA = itr
        clrA = _zone_color(A)
        mkrA = _zone_markersize(A)
        scatter!(ax, r[itrA], X[itrA], markersize = mkrA, color=clrA)
    elseif B == A + 1
        itrA = A == 1 ? (N1:k+1) : A == 2 ? (N1:Na) : A == 3 ? (N1:Nuctp) : A == 4 ? (N1:Nb) : A == 5 ? (N1:N-k-1) : 0:0
        itrB = B == 2 ? (k+2:N2) : B == 3 ? (Na+1:N2) : B == 4 ? (Nuctp+1:N2) : B == 5 ? (Nb+1:N2) : B == 6 ? (N-k:N2) : 0:0
        clrA = _zone_color(A)
        clrB = _zone_color(B)
        mkrA = _zone_markersize(A)
        mkrB = _zone_markersize(B)
        scatter!(ax, r[itrA], X[itrA], markersize = mkrA, color=clrA)
        scatter!(ax, r[itrB], X[itrB], markersize = mkrB, color=clrB)
    elseif B == A + 2
        C = A + 1
        itrA = A == 1 ? (N1:k+1)  : A == 2 ? (N1:Na)      : A == 3 ? (N1:Nuctp)   : A == 4 ? (N1:Nb)      : 0:0
        itrC = C == 2 ? (k+2:Na)  : C == 3 ? (Na+1:Nuctp) : C == 4 ? (Nuctp+1:Nb) : C == 5 ? (Nb+1:N-k-1) : 0:0
        itrB = B == 3 ? (Na+1:N2) : B == 4 ? (Nuctp+1:N2) : B == 5 ? (Nb+1:N2)    : B == 6 ? (N-k:N2)     : 0:0
        clrA = _zone_color(A)
        clrB = _zone_color(B)
        clrC = _zone_color(C)
        mkrA = _zone_markersize(A)
        mkrB = _zone_markersize(B)
        mkrC = _zone_markersize(C)
        scatter!(ax, r[itrA], X[itrA], markersize = mkrA, color=clrA)
        scatter!(ax, r[itrB], X[itrB], markersize = mkrB, color=clrB)
        scatter!(ax, r[itrC], X[itrC], markersize = mkrC, color=clrC)
    elseif B == A + 3
        C = A + 1
        D = A + 2
        itrA = A == 1 ? (N1:k+1)     : A == 2 ? (N1:Na)      : A == 3 ? (N1:Nuctp)   : 0:0
        itrC = C == 2 ? (k+2:Na)     : C == 3 ? (Na+1:Nuctp) : C == 4 ? (Nuctp+1:Nb) : 0:0
        itrD = D == 3 ? (Na+1:Nuctp) : D == 4 ? (Nuctp+1:Nb) : D == 5 ? (Nb+1:N-k-1) : 0:0
        itrB = B == 4 ? (Nuctp+1:N2) : B == 5 ? (Nb+1:N2) : B == 6 ? (N-k:N2) : 0:0
        clrA = _zone_color(A)
        clrB = _zone_color(B)
        clrC = _zone_color(C)
        clrD = _zone_color(D)
        mkrA = _zone_markersize(A)
        mkrB = _zone_markersize(B)
        mkrC = _zone_markersize(C)
        mkrD = _zone_markersize(D)
        scatter!(ax, r[itrA], X[itrA], markersize = mkrA, color=clrA)
        scatter!(ax, r[itrB], X[itrB], markersize = mkrB, color=clrB)
        scatter!(ax, r[itrC], X[itrC], markersize = mkrC, color=clrC)
        scatter!(ax, r[itrD], X[itrD], markersize = mkrD, color=clrD)
    elseif B == A + 4
        C = A + 1
        D = A + 2
        E = A + 3
        itrA = A == 1 ? (N1:k+1)     : A == 2 ? (N1:Na)      : 0:0
        itrC = C == 2 ? (k+2:Na)     : C == 3 ? (Na+1:Nuctp) : 0:0
        itrD = D == 3 ? (Na+1:Nuctp) : D == 4 ? (Nuctp+1:Nb) : 0:0
        itrE = E == 4 ? (Nuctp+1:Nb) : E == 5 ? (Nb+1:N-k-1) : 0:0
        itrB = B == 5 ? (Nb+1:N2)    : B == 6 ? (N-k:N2)     : 0:0
        clrA = _zone_color(A)
        clrB = _zone_color(B)
        clrC = _zone_color(C)
        clrD = _zone_color(D)
        clrE = _zone_color(E)
        mkrA = _zone_markersize(A)
        mkrB = _zone_markersize(B)
        mkrC = _zone_markersize(C)
        mkrD = _zone_markersize(D)
        mkrE = _zone_markersize(E)
        scatter!(ax, r[itrA], X[itrA], markersize = mkrA, color=clrA)
        scatter!(ax, r[itrB], X[itrB], markersize = mkrB, color=clrB)
        scatter!(ax, r[itrC], X[itrC], markersize = mkrC, color=clrC)
        scatter!(ax, r[itrD], X[itrD], markersize = mkrD, color=clrD)
        scatter!(ax, r[itrE], X[itrE], markersize = mkrE, color=clrE)
    else #    (B = A + 5)
        C = A + 1
        D = A + 2
        E = A + 3
        F = A + 4
        itrA = A == 1 ? (N1:k+1)     : 0:0
        itrC = C == 2 ? (k+2:Na)     : 0:0
        itrD = D == 3 ? (Na+1:Nuctp) : 0:0
        itrE = E == 4 ? (Nuctp+1:Nb) : 0:0
        itrF = F == 5 ? (Nb+1:N-k-1) : 0:0
        itrB = B == 6 ? (N-k:N2)     : 0:0
        clrA = _zone_color(A)
        clrB = _zone_color(B)
        clrC = _zone_color(C)
        clrD = _zone_color(D)
        clrE = _zone_color(E)
        clrF = _zone_color(F)
        mkrA = _zone_markersize(A)
        mkrB = _zone_markersize(B)
        mkrC = _zone_markersize(C)
        mkrD = _zone_markersize(D)
        mkrE = _zone_markersize(E)
        mkrF = _zone_markersize(F)
        scatter!(ax, r[itrA], X[itrA], markersize = mkrA, color=clrA)
        scatter!(ax, r[itrB], X[itrB], markersize = mkrB, color=clrB)
        scatter!(ax, r[itrC], X[itrC], markersize = mkrC, color=clrC)
        scatter!(ax, r[itrD], X[itrD], markersize = mkrD, color=clrD)
        scatter!(ax, r[itrE], X[itrE], markersize = mkrE, color=clrE)
        scatter!(ax, r[itrF], X[itrF], markersize = mkrF, color=clrF)
    end

    N1 ≤ Nuctp ≤ N2 ? scatter!(ax, (r[Nuctp], X[Nuctp]), markersize = 5, color=:Black) : false

end

# =========================================================================================================================

function plot_OUTSCH(Z::Vector{Complex{T}}, grid::Grid{T}, def::Def{T}; theme = Theme(fontsize = 10, resolution = (900,900)) ) where T<:Real

    set_theme!(theme)    
    
    r = grid.r
    k = def.k
    l = k+1

       Na = def.pos.Na
     Nmin = def.pos.Nmin
    Nuctp = def.pos.Nuctp
        N = def.pos.N
       Nb = def.pos.Nb

    itr = 1:k+1

    (ylabelP1, ylabelQ1, Z1) = ("χ(r)", "dχ/dr", Z)
    (ylabelP2, ylabelQ2, Z2) = ("ψ(r)", "dψ/dr", wavefunction(Z, grid))


    R = [r[n] for n ∈ itr]
    P1 = [real(Z1[n]) for n ∈ itr]
    Q1 = [imag(Z1[n]) for n ∈ itr]
    P2 = [real(Z2[n]) for n ∈ itr]
    Q2 = [imag(Z2[n]) for n ∈ itr]

    symbol = def.atom.element.symbol
      name = def.orbit.name

    fig = Figure()

    attr1a = set_attributes(fig; title = "OUTSCH: " * ylabelP1 * " on Grid[$(itr)]", xlabel = "r (a.u.)", ylabel = ylabelP1 )
    attr1b = set_attributes(fig; title = "OUTSCH: " * ylabelQ1 * " on Grid[$(itr)]", xlabel = "r (a.u.)", ylabel = ylabelQ1 )
    attr2a = set_attributes(fig; title = "OUTSCH: " * ylabelP2 * " on Grid[$(itr)]", xlabel = "r (a.u.)", ylabel = ylabelP2 )
    attr2b = set_attributes(fig; title = "OUTSCH: " * ylabelQ2 * " on Grid[$(itr)]", xlabel = "r (a.u.)", ylabel = ylabelQ2 )

    ax0 = Label(fig; text = symbol * ": " * name, textsize = 24,  color=:gray)
    ax1a = Axis(fig; attr1a..., yticks = P1)                   # create axes, add attributes
    ax1b = Axis(fig; attr1b..., yticks = Q1)                   # create axes, add attributes 
    ax2a = Axis(fig; attr2a..., yticks = P2)                   # create axes, add attributes
    ax2b = Axis(fig; attr2b..., yticks = Q2)                   # create axes, add attributes 

    lines!(ax1a, R  , P1, markersize = 1, color=:gray90)
    scatter!(ax1a, R, P1, markersize = 2, color=:black)
    lines!(ax1b, R  , Q1, markersize = 1, color=:gray90)
    scatter!(ax1b, R, Q1, markersize = 2, color=:black)
    lines!(ax2a, R  , P2, markersize = 1, color=:gray90)
    scatter!(ax2a, R, P2, markersize = 2, color=:black)
    lines!(ax2b, R  , Q2, markersize = 1, color=:gray90)
    scatter!(ax2b, R, Q2, markersize = 2, color=:black)

    fig[1,1] = ax1a                                      # create layout and show figure
    fig[1,2] = ax1b                                      # create layout and show figure
    fig[2,1] = ax2a                                      # create layout and show figure
    fig[2,2] = ax2b                                      # create layout and show figure
    fig[0,:] = ax0

    return fig

end

# =========================================================================================================================

function plot_endpoints(Z::Vector{Complex{T}}, grid::Grid{T}, def::Def{T}; resolution = (900,600)) where T<:Real

    theme = Theme(fontsize = 10, resolution)
    set_theme!(theme)

    r = grid.r
    k = def.k

       Na = def.pos.Na
     Nmin = def.pos.Nmin
    Nuctp = def.pos.Nuctp
        N = def.pos.N
       Nb = def.pos.Nb

    itr1 = 1:Na+k
    itr2 = Nb-k:N

    symbol = def.atom.element.symbol
      name = def.orbit.name

    fig = Figure()

    attr1 = set_attributes(fig; title = "OUTSCH: χ(r) on Grid[$(itr1)]", xlabel = "r (a.u.)", ylabel = "χ(r)")
    attr2 = set_attributes(fig; title = "INSCH:  χ(r) on Grid[$(itr2)]", xlabel = "r (a.u.)", ylabel = "χ(r)")
    attr3 = set_attributes(fig; title = "OUTSCH: dχ/dr on Grid[$(itr1)]", xlabel = "r (a.u.)", ylabel = "dχ/dr")
    attr4 = set_attributes(fig; title = "INSCH:  dχ/dr on Grid[$(itr2)]", xlabel = "r (a.u.)", ylabel = "dχ/dr")

    ax0 = Label(fig; text = symbol * ": " * name * " - Endpoints", textsize = 24,  color=:gray)
    ax1 = Axis(fig; attr1...)                            # create axes, add atrributes
    ax2 = Axis(fig; attr2...)                            # create axes, add atrributes
    ax3 = Axis(fig; attr3...)                            # create axes, add atrributes
    ax4 = Axis(fig; attr4...)                            # create axes, add atrributes

      lines!(ax1, r[itr1], real(Z[itr1]), markersize = 1, color=:gray90)
    _scatterX!(ax1, r, real(Z), itr1, def)
      lines!(ax2, r[itr2], real(Z[itr2]), markersize = 1, color=:gray90)
    _scatterX!(ax2, r, real(Z), itr2, def)
      lines!(ax3, r[itr1], imag(Z[itr1]), markersize = 2, color=:gray90)
    _scatterX!(ax3, r, imag(Z), itr1, def)
      lines!(ax4, r[itr2], imag(Z[itr2]), markersize = 2, color=:gray90)
    _scatterX!(ax4, r, imag(Z), itr2, def)

    fig[1,1] = ax1                                      # create layout and show figure
    fig[1,2] = ax2                                      # create layout and show figure
    fig[2,1] = ax3                                      # create layout and show figure
    fig[2,2] = ax4                                      # create layout and show figure
    fig[0,:] = ax0

    return fig

end

# ===============================================================================
function plot_wavefunction(Nstart::Int, Nstop::Int, E::T, grid::Grid{T}, def::Def{T}, Z::Vector{Complex{T}}; theme = Theme(fontsize = 10, resolution = (900,600)), reduced=true) where T <: Real

    set_theme!(theme)

        r = grid.r
        k = def.k
       Na = def.pos.Na
    Nuctp = def.pos.Nuctp
       Nb = def.pos.Nb
        N = def.pos.N

    println()
    println("plot_wavefunction: color coding based on k+1 = $(k+1), Na = $(def.pos.Na), Nuctp = $(def.pos.Nuctp), Nb = $(def.pos.Nb), N-k = $(def.pos.N-k), N = $(def.pos.N)")

    1 ≤ Nstart ≤ N || error("Error: Nstart outside index range 1:$(N)")
    1 ≤ Nstop  ≤ N || error("Error: Nstop outside index range 1:$(N)")

    ylabel1, ylabel2, Z = reduced ? ("χ(r)", "dχ/dr", Z) :  ("ψ(r)", "dψ/dr", wavefunction(Z, grid))
    ylabelA, ylabelB = reduced ? ("χ(n)", "dχ/dr") :  ("ψ(n)", "dψ/dr")

    header = reduced ? "reduced wavefunction" :  "full wavefunction"

    A = _zone(Nstart, def)
    B = _zone(Nstop, def)

    P = convert.(Float64,real(Z))
    Q = convert.(Float64,imag(Z))

    P0 = P
    Q0 = Q

    r = convert.(Float64,r)
    n = [i for i=1:N]

    itr0 = Nstart:Nstop

    symbol = def.atom.element.symbol
      name = def.orbit.name

    fig = Figure()

    attr1 = set_attributes(fig; title = ylabel1 * " on Grid[$(itr0)]", xlabel = "r (a.u.)", ylabel = ylabel1)
    attr2 = set_attributes(fig; title = ylabelA * " on Grid[$(itr0)]", xlabel = "n", ylabel = ylabelA)
    attr3 = set_attributes(fig; title = ylabel2 * " on Grid[$(itr0)]", xlabel = "r (a.u.)", ylabel = ylabel2)
    attr4 = set_attributes(fig; title = ylabelB * " on Grid[$(itr0)]", xlabel = "n", ylabel = ylabelB)

    ax0 = Label(fig; text = symbol * ": " * name, textsize = 24,  color=:gray)
    ax1 = Axis(fig; attr1...)
    ax2 = Axis(fig; attr2...)
    ax3 = Axis(fig; attr3...)
    ax4 = Axis(fig; attr4...)

    lines!(ax1, r[itr0], P[itr0] , markersize = 1, color=:gray90)
    lines!(ax2, n[itr0], P0[itr0], markersize = 1, color=:gray90)
    lines!(ax3, r[itr0], Q[itr0] , markersize = 1, color=:gray90)
    lines!(ax4, n[itr0], Q0[itr0], markersize = 1, color=:gray90)

    lines!(ax1, r[itr0], P[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax1, r, P, itr0, def)

    lines!(ax2, n[itr0], P0[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax2, n, P0, itr0, def)

    lines!(ax3, r[itr0], Q[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax3, r, Q, itr0, def)

    lines!(ax4, n[itr0], Q0[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax4, n, Q0, itr0, def)

    fig[1,1] = ax1
    fig[1,2] = ax2
    fig[2,1] = ax3
    fig[2,2] = ax4
    fig[0,:] = ax0

    return fig

end
function plot_wavefunction(itr::UnitRange{Int}, E::T, grid::Grid{T}, def::Def{T}, Z::Vector{Complex{T}}; theme = Theme(fontsize = 10, resolution = (900,600)), reduced=true) where T <: Real

    Nstart = itr.start
    Nstop  = itr.stop

    return plot_wavefunction(Nstart, Nstop, E, grid, def, Z; theme, reduced)

end

# ===============================================================================
function compare_wavefunctions(Z1::Vector{Complex{T}}, Z2::Vector{Complex{T}}, Nstart::Int, Nstop::Int, E::T, grid::Grid{T}, def::Def{T}; theme = Theme(fontsize = 10, resolution = (900,900)), reduced=true) where T <: Real

    set_theme!(theme)

        r = grid.r
        k = def.k
       Na = def.pos.Na
    Nuctp = def.pos.Nuctp
       Nb = def.pos.Nb
        N = def.pos.N

    println()
    println("wavefunction color coding based on k+1 = $(k+1), Na = $(def.pos.Na), Nuctp = $(def.pos.Nuctp), Nb = $(def.pos.Nb), N-k = $(def.pos.N-k), N = $(def.pos.N)")

    1 ≤ Nstart ≤ N || error("Error: Nstart outside index range 1:$(N)")
    1 ≤ Nstop  ≤ N || error("Error: Nstop outside index range 1:$(N)")

    ylabel1a, ylabel1b, Z1 = reduced ? ("χ1(r)", "dχ1/dr", Z1) :  ("ψ1(r)", "dψ1/dr", wavefunction(Z1, grid))
    ylabel2a, ylabel2b, Z2 = reduced ? ("χ2(r)", "dχ2/dr", Z2) :  ("ψ2(r)", "dψ2/dr", wavefunction(Z2, grid))
    ylabel3a, ylabel3b, Z3 = reduced ? ("χ2(r)-χ1(r)", "d(χ2-χ1)/dr", Z2) :  ("ψ2(r)-ψ1(r)", "d(ψ2-ψ1)/dr", wavefunction(Z3, grid))
    ylabel4a, ylabel4b, Z4 = reduced ? ("χ2(r)/χ1(r)-1", "d(χ2/χ1)/dr", Z2) :  ("ψ2(r)/ψ1(r)", "d(ψ2/ψ1)/dr", wavefunction(Z3, grid))
    ylabel1A, ylabel1B = reduced ? ("χ1(n)", "dχ1/dr") :  ("ψ1(n)", "dψ1/dr")
    ylabel2A, ylabel2B = reduced ? ("χ2(n)", "dχ2/dr") :  ("ψ2(n)", "dψ2/dr")
    ylabel3A, ylabel3B = reduced ? ("χ2(n)-χ1(r)", "d(χ2-χ1)/dr") :  ("ψ2(n)-ψ1(n)", "d(ψ2-ψ1)/dr")
    ylabel4A, ylabel4B = reduced ? ("χ2(n)/χ1(r)-1", "d(χ2/χ1)/dr") :  ("ψ2(n)/ψ1(n)-1", "d(ψ2/ψ1)/dr")

    header = reduced ? "reduced wavefunction" :  "full wavefunction"

    A = _zone(Nstart, def)
    B = _zone(Nstop, def)

    P1 = convert.(Float64,real(Z1))
    Q1 = convert.(Float64,imag(Z1))
    
    P2 = convert.(Float64,real(Z2))
    Q2 = convert.(Float64,imag(Z2))
    
    P3 = P2 .- P1
    Q3 = Q2 .- Q1
    
    P4 = P2 ./ P1 .- 1
    Q4 = Q2 ./ Q1 .- 1
    
    P4[1] = ((P1[1] == 0) | (P2[1] == 0)) ? fdiff_interpolation(P4[2:k+1], 0) : P4[1]
    Q4[1] = ((Q1[1] == 0) | (Q2[1] == 0)) ? fdiff_interpolation(Q4[2:k+1], 0) : Q4[1]
    
    P10 = P1
    Q10 = Q1
    
    P20 = P2
    Q20 = Q2
    
    P30 = P3
    Q30 = Q3
    
    P40 = P4
    Q40 = Q4
    

    r = convert.(Float64,r)
    n = [i for i=1:N]

    itr0 = Nstart:Nstop

    symbol = def.atom.element.symbol
      name = def.orbit.name

    fig = Figure()

    attr1a = set_attributes(fig; title = ylabel1a, xlabel = "r (a.u.)", ylabel = ylabel1a)
    attr1A = set_attributes(fig; title = ylabel1A, xlabel = "n", ylabel = ylabel1A)
    attr1b = set_attributes(fig; title = ylabel1b, xlabel = "r (a.u.)", ylabel = ylabel1b)
    attr1B = set_attributes(fig; title = ylabel1B, xlabel = "n", ylabel = ylabel1B)
    attr2a = set_attributes(fig; title = ylabel2a, xlabel = "r (a.u.)", ylabel = ylabel2a)
    attr2A = set_attributes(fig; title = ylabel2A, xlabel = "n", ylabel = ylabel2A)
    attr2b = set_attributes(fig; title = ylabel2b, xlabel = "r (a.u.)", ylabel = ylabel2b)
    attr2B = set_attributes(fig; title = ylabel2B, xlabel = "n", ylabel = ylabel2B)
    attr3a = set_attributes(fig; title = ylabel3a, xlabel = "r (a.u.)", ylabel = ylabel3a)
    attr3A = set_attributes(fig; title = ylabel3A, xlabel = "n", ylabel = ylabel3A)
    attr3b = set_attributes(fig; title = ylabel3b, xlabel = "r (a.u.)", ylabel = ylabel3b)
    attr3B = set_attributes(fig; title = ylabel3B, xlabel = "n", ylabel = ylabel3B)
    attr4a = set_attributes(fig; title = ylabel4a, xlabel = "r (a.u.)", ylabel = ylabel4a)
    attr4A = set_attributes(fig; title = ylabel4A, xlabel = "n", ylabel = ylabel4A)
    attr4b = set_attributes(fig; title = ylabel4b, xlabel = "r (a.u.)", ylabel = ylabel4b)
    attr4B = set_attributes(fig; title = ylabel4B, xlabel = "n", ylabel = ylabel4B)

    ax0 = Label(fig; text = symbol * ": " * name * " on Grid[$(itr0)]", textsize = 24,  color=:gray)
    ax1a = Axis(fig; attr1a...)
    ax1A = Axis(fig; attr1A...)
    ax1b = Axis(fig; attr1b...)
    ax1B = Axis(fig; attr1B...)
    ax2a = Axis(fig; attr2a...)
    ax2A = Axis(fig; attr2A...)
    ax2b = Axis(fig; attr2b...)
    ax2B = Axis(fig; attr2B...)
    ax3a = Axis(fig; attr3a...)
    ax3A = Axis(fig; attr3A...)
    ax3b = Axis(fig; attr3b...)
    ax3B = Axis(fig; attr3B...)
    ax4a = Axis(fig; attr4a...)
    ax4A = Axis(fig; attr4A...)
    ax4b = Axis(fig; attr4b...)
    ax4B = Axis(fig; attr4B...)


    lines!(ax1a, r[itr0], P1[itr0] , markersize = 1, color=:gray90)
    lines!(ax1A, n[itr0], P10[itr0], markersize = 1, color=:gray90)
    lines!(ax1b, r[itr0], Q1[itr0] , markersize = 1, color=:gray90)
    lines!(ax1B, n[itr0], Q10[itr0], markersize = 1, color=:gray90)
    lines!(ax2a, r[itr0], P2[itr0] , markersize = 1, color=:gray90)
    lines!(ax2A, n[itr0], P20[itr0], markersize = 1, color=:gray90)
    lines!(ax2b, r[itr0], Q2[itr0] , markersize = 1, color=:gray90)
    lines!(ax2B, n[itr0], Q20[itr0], markersize = 1, color=:gray90)
    lines!(ax3a, r[itr0], P3[itr0] , markersize = 1, color=:gray90)
    lines!(ax3A, n[itr0], P30[itr0], markersize = 1, color=:gray90)
    lines!(ax3b, r[itr0], Q3[itr0] , markersize = 1, color=:gray90)
    lines!(ax3B, n[itr0], Q30[itr0], markersize = 1, color=:gray90)
    lines!(ax4a, r[itr0], P4[itr0] , markersize = 1, color=:gray90)
    lines!(ax4A, n[itr0], P40[itr0], markersize = 1, color=:gray90)
    lines!(ax4b, r[itr0], Q4[itr0] , markersize = 1, color=:gray90)
    lines!(ax4B, n[itr0], Q40[itr0], markersize = 1, color=:gray90)

    lines!(ax1a, r[itr0], P1[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax1a, r, P1, itr0, def)

    lines!(ax1A, n[itr0], P10[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax1A, n, P10, itr0, def)

    lines!(ax1b, r[itr0], Q1[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax1b, r, Q1, itr0, def)

    lines!(ax1B, n[itr0], Q10[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax1B, n, Q10, itr0, def)
    
    lines!(ax2a, r[itr0], P2[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax2a, r, P2, itr0, def)

    lines!(ax2A, n[itr0], P20[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax2A, n, P20, itr0, def)

    lines!(ax2b, r[itr0], Q2[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax2b, r, Q2, itr0, def)

    lines!(ax2B, n[itr0], Q20[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax2B, n, Q20, itr0, def)
    
    lines!(ax3a, r[itr0], P3[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax3a, r, P3, itr0, def)

    lines!(ax3A, n[itr0], P30[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax3A, n, P30, itr0, def)

    lines!(ax3b, r[itr0], Q3[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax3b, r, Q3, itr0, def)

    lines!(ax3B, n[itr0], Q30[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax3B, n, Q30, itr0, def)
    
    lines!(ax4a, r[itr0], P4[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax4a, r, P4, itr0, def)

    lines!(ax4A, n[itr0], P40[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax4A, n, P40, itr0, def)

    lines!(ax4b, r[itr0], Q4[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax4b, r, Q4, itr0, def)

    lines!(ax4B, n[itr0], Q40[itr0], markersize = 1, color=:gray90)
    _scatterX!(ax4B, n, Q40, itr0, def)

    fig[1,1] = ax1a
    fig[2,1] = ax2a
    fig[3,1] = ax3a
    fig[4,1] = ax4a
    
    fig[1,2] = ax1A
    fig[2,2] = ax2A
    fig[3,2] = ax3A
    fig[4,2] = ax4A
 
    fig[1,3] = ax1b
    fig[2,3] = ax2b   
    fig[3,3] = ax3b   
    fig[4,3] = ax4b
    
    fig[1,4] = ax1B
    fig[2,4] = ax2B
    fig[3,4] = ax3B
    fig[4,4] = ax4B
    
    fig[0,:] = ax0

    return fig

end
function compare_wavefunctions(Z1::Vector{Complex{T}}, Z2::Vector{Complex{T}}, itr::UnitRange{Int}, E::T, grid::Grid{T}, def::Def{T}; theme = Theme(fontsize = 10, resolution = (900,900)), reduced=true) where T <: Real

    Nstart = itr.start
    Nstop  = itr.stop

    return compare_wavefunctions(Z1, Z2, Nstart, Nstop, E, grid, def; theme, reduced)

end

# =========================================================================================================================

function analyze_wavefunction(Z::Vector{Complex{T}}, grid::Grid{T}, def::Def{T}; theme = Theme(fontsize = 10, resolution = (900,600)), reduced=true) where T <: Real

    set_theme!(theme)

        r = grid.r
      pos = def.pos
       Na = def.pos.Na
    Nuctp = def.pos.Nuctp
       Nb = def.pos.Nb
        N = def.pos.N

    ylabel1, ylabel2, Z = reduced ? ("χ(r)", "dχ/dr", Z) :  ("ψ(r)", "dψ/dr", wavefunction(Z,grid))

    header = reduced ? "reduced wavefunction" :  "full wavefunction"
    P = real(Z)
    Q = imag(Z)

    itr0 = 1:N
    itr1 = 1:Nuctp
    itr2 = Nuctp-9:Nuctp+9
    itr3 = Nuctp+1:N
    itr4 = 1:Na+10
    itr5 = itr2
    itr6 = Nb-9:N

    symbol = def.atom.element.symbol
      name = def.orbit.name

    fig = Figure()

    attr1 = set_attributes(fig; title = "Grid[$(itr0)]", xlabel = "r (a.u.)", ylabel = ylabel1)
    attr2 = set_attributes(fig; title = "Grid[$(itr2)]", xlabel = "r (a.u.)", ylabel = ylabel1)
    attr3 = set_attributes(fig; title = "Grid[$(itr0)]", xlabel = "r (a.u.)", ylabel = ylabel2)
    attr4 = set_attributes(fig; title = "Grid[$(itr4)]", xlabel = "r (a.u.)", ylabel = ylabel1)
    attr5 = set_attributes(fig; title = "Grid[$(itr5)]", xlabel = "r (a.u.)", ylabel = ylabel2)
    attr6 = set_attributes(fig; title = "Grid[$(itr6)]", xlabel = "r (a.u.)", ylabel = ylabel1)

    ax0 = Label(fig; text = symbol * ": " * name * " - " * header, textsize = 24,  color=:gray)
    ax1 = Axis(fig; attr1...)
    ax2 = Axis(fig; attr2...)
    ax3 = Axis(fig; attr3...)
    ax4 = Axis(fig; attr4...)
    ax5 = Axis(fig; attr5...)
    ax6 = Axis(fig; attr6...)

    lines!(ax1, r, P, markersize = 1, color=:gray90)
    _scatterX!(ax1, r, P, itr0, def)

    lines!(ax2, r[itr2], P[itr2], markersize = 1, color=:gray90)
    _scatterX!(ax2, r, P, itr2, def)

    lines!(ax3, r, Q, markersize = 1, color=:gray90)
    _scatterX!(ax3, r, Q, itr0, def)

    lines!(ax4, r[itr4], P[itr4]/P[Na], markersize = 1, color=:gray90)
    _scatterX!(ax4, r, P/P[Na], itr4, def)

    lines!(ax5, r[itr5], Q[itr5], markersize = 1, color=:gray90)
    _scatterX!(ax5, r, Q, itr5, def)

    lines!(ax6, r[itr6], P[itr6], markersize = 1, color=:gray90)
    _scatterX!(ax6, r, P, itr6, def)

    fig[1,1] = ax1
    fig[1,2] = ax2
    fig[1,3] = ax3
    fig[2,1] = ax4
    fig[2,2] = ax5
    fig[2,3] = ax6
    fig[0,:] = ax0

    return fig

end

println("Included:
plot_tabulated_function(Nstart, Nstop, f, grid; title, theme)
plot_gridfunction(Nstart, Nstop, grid; title='')
plot_potentials(E, grid, def)
plot_wavefunction(Nstart, Nstop, E, grid, def, Z; reduced=true)
compare_wavefunctions(Z1, Z2, Nstart, Nstop, E, grid, def; reduced=true)
plot_endpoints(Z, grid, def)
analyze_wavefunction(Z, grid, def; reduced=true)")
