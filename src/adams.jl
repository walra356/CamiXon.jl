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
   castAdams(E::T, grid::Grid{T}, def::Def{T}) where T<:Real

Initiates [`Adams`](@ref) object.
"""
function castAdams(E::T, grid::Grid{T}, def::Def{T}) where T<:Real

    castPos(E, def.potscr, grid)

    G = matG(E, grid, def)
    σ = matσ(E, grid, def)
    M = matMinv(E, grid, def)

    return Adams1(G, σ, M)

end

# ------------------------------------------------------------------------------
#                   updateAdams!(adams, E, grid, def)
# ------------------------------------------------------------------------------


@doc raw"""
    updateAdams!(adams::Adams1{T}, E::T, grid::Grid{T}, def::Def{T}) where T<:Real

Upate [`Adams`](@ref) object.
"""
function updateAdams!(adams::Adams1{T}, E::T, grid::Grid{T}, def::Def{T}) where T<:Real

G = matG(E, grid, def)
σ = matσ(E, grid, def)
M = matMinv(E, grid, def)

return Adams1(G, σ, M)

end
