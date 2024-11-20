# SPDX-License-Identifier: MIT

# author: Jook Walraven - 14-2-2023

# ==============================================================================
#                               insch.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#                       INSCH!(Z, E, grid, def)
# ------------------------------------------------------------------------------

@doc raw"""
    INSCH!(Z::Vector{Complex{T}}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real
 
Ansatz solution for the *inward* integration of the radial wave equation for the first ``k`` points 
on the [`CamiDiff.Grid`](@ref), where ``k`` is the Adams-Moulton order. The Ansatz is based on the WKB solution 
for energy `E` at distances *far above* the upper classical turning point - uctp)
"""
function INSCH!(Z::Vector{Complex{T}}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real
    
    Z = INSCH_WKB!(Z, E, grid, def)

    return Z

end

# ------------------------------------------------------------------------------
#                       INSCH_WKB!(Z, E, grid, def)
# ------------------------------------------------------------------------------

@doc raw"""
    INSCH_WKB!(Z::Vector{Complex{T}}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real

WKB Ansatz of `k+1` points for [`INSCH!`](@ref)
"""
function INSCH_WKB!(Z::Vector{Complex{T}}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real

    N = grid.N
    r = grid.r
    k = grid.k
    Nuctp = def.pos.Nuctp
    two = T(2)
    pot = def.potscr

    p = Array{T,1}(undef,N)
    I = Array{T,1}(undef,N)
    P = Array{T,1}(undef,N)
    Q = Array{T,1}(undef,N)
    
    for n=Nuctp:N
        p[n] = sqrt(abs(pot[n]-E))                         # quasi-classical momentum   
        I[n] = CamiDiff.grid_integration(p, grid, Nuctp:n)
        P[n] = exp(-I[n])/sqrt(p[n])                       # WKB solution
        def.pos.Nb = P[n] > 1.0e-30 ? n : break
    end

    Nb = def.pos.Nb = def.pos.Nb - k
    
    Q[Nb-k:Nb+k] = CamiDiff.grid_differentiation(P, grid, Nb-k:Nb+k)   # avoid lower end point correction by doubling range
    
    for n=Nb-k:Nb+k
        Z[n] = P[n] + im * Q[n]
    end 
    
    return Z

end

# ------------------------------------------------------------------------------
#                       INSCH_WJ!(Z, E, grid, def) 
# kept for the record in "assets" folder
# ------------------------------------------------------------------------------