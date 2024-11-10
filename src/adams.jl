# SPDX-License-Identifier: MIT

# author: Jook Walraven - 4-9-2024

# ==============================================================================
#                               adams.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#                           Adams(G, σ, Minv) 
# ------------------------------------------------------------------------------

@doc raw"""
    Adams{T}

* G: (`:Vector{Matrix{T}}`)
* σ: (`:Vector{Matrix{T}}`)
* Minv: (`:Vector{Matrix{T}}`)
"""
struct Adams{T}

    G::Vector{Matrix{T}}
    σ::Vector{Matrix{T}}
    Minv::Vector{Matrix{T}}
    Z::Vector{Complex{T}} 

end
struct Adams1{T}

    G::Vector{Matrix{T}}
    σ::Vector{Matrix{T}}
    Minv::Vector{Matrix{T}}

end

# ------------------------------------------------------------------------------
#                   castAdams(E, grid, def)
# ------------------------------------------------------------------------------

@doc raw"""
    castAdams1(E::T, grid::Grid{T}, def::Def{T}) where T<:Real

Initiates [`Adams`](@ref) object.
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
    Z = OUTSCH(E, grid, def, σ)

    return Adams(G, σ, M, Z)

end
function castAdams1(E::T, grid::Grid{T}, def::Def{T}) where T<:Real

    castPos(E, def.potscr, grid)

    G = matG1(E, grid, def)
    σ = matσ1(E, grid, def)
    M = matMinv1(E, grid, def)

    return Adams1(G, σ, M)

end

# ------------------------------------------------------------------------------
#                   updateAdams!(adams, E, grid, def)
# ------------------------------------------------------------------------------


@doc raw"""
    updateAdams(adams::Adams{T}, E::T, grid::Grid{T}, def::Def{T}) where T<:Real

Upate [`Adams`](@ref) object.
"""
function updateAdams!(adams::Adams{T}, E::T, grid::Grid{T}, def::Def{T}) where T<:Real

    E = convert(T, E)
    G = matG(E, grid, def)
    σ = matσ(E, grid, def)
    M = matMinv(E, grid, def, def.am[end])
    Z = adams.Z

    def.pos.Na = get_Na(Z, def)
    def.pos.Nlctp = get_Nlctp(E, def)
    def.pos.Nuctp = get_Nuctp(E, def)

    return Adams(G, σ, M, Z)

end
function updateAdams1(adams::Adams1{T}, E::T, grid::Grid{T}, def::Def{T}) where T<:Real

G = matG1(E, grid, def)
σ = matσ1(E, grid, def)
M = matMinv1(E, grid, def)

return Adams1(G, σ, M)

end
