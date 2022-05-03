# =========================== Adams sector =====================================

@doc raw"""
    Adams

* G: (`:Vector{Matrix{T}}`)
* σ: (`:Vector{Matrix{T}}`)
* Minv: (`:Vector{Matrix{T}}`)
* Z: (`:Vector{Complex{T}}`)
"""
struct Adams{T}

    G::Vector{Matrix{T}}
    σ::Vector{Matrix{T}}
    Minv::Vector{Matrix{T}}
    Z::Vector{Complex{T}}

end

@doc raw"""
    castAdams(E::T, grid::Grid{T}, def::Def{T}) where T<:Real

"""
function castAdams(E::T, grid::Grid{T}, def::Def{T}) where T<:Real

    am = def.am
     k = def.k
     ℓ = def.orbit.ℓ

    N = def.pos.N
    def.pos.Nmin = get_Nmin(def)
    def.pos.Nlctp = get_Nlctp(E, def)
    def.pos.Nuctp = get_Nuctp(E, def)

    G = matG(E, grid, def)
    σ = matσ(E, grid, def)
    M = matMinv(E, grid, def, am[end])
    Z = ℓ < 5 ? OUTSCH(grid, def, σ) : OUTSCH_WKB(E, grid, def)

    return Adams(G, σ, M, Z)

end

@doc raw"""
    updateAdams!(adams::Adams{T}, E, grid::Grid{T}, def::Def{T}) where T<:Real

"""
function updateAdams!(adams::Adams{T}, E, grid::Grid{T}, def::Def{T}) where T<:Real

    E = convert(def.T, E)
    G = matG(E, grid, def)
    σ = matσ(E, grid, def)
    M = matMinv(E, grid, def, def.am[end])
    Z = adams.Z

    def.pos.Na = get_Na(Z, def)
    def.pos.Nlctp = get_Nlctp(E, def)
    def.pos.Nuctp = get_Nuctp(E, def)

    return Adams(G, σ, M, Z)

end
