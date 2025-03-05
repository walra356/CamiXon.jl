# SPDX-License-Identifier: MIT

# author: Jook Walraven - 11-11-2024

# ==============================================================================
#                               init.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#                               Init{T}
# ------------------------------------------------------------------------------

@doc raw"""
    Init{T} where T<:Real

Energy object with fields:

* ` .Emin` : lower energy limit (`::T`)
* ` .E` : trial energy (`::T`)
* ` .Emax` : upper energy limit (`::T`)
* ` .ΔE` : difference with respect to previous E used to calculate current E  (`::T`)
"""
mutable struct Init{T}

    Emin::T     # lowest allowed energy
    E::T        # trial energy
    Emax::T     # practical highest energy
    ΔE::T       # ΔE with respect to previous E to calculate E trial in use
    
end

# ------------------------------------------------------------------------------
#                               castInit{T}
# ------------------------------------------------------------------------------

@doc raw"""
    castInit(E::T, def::Def{T}) where T<:Real

Create and initialize the energy object [`Init`](@ref) according to the 
specifications given in the [`Def`](@ref) object.

* ` .Emin` : lower energy limit (`::T`)
* ` .E` : trial energy (`::T`)
* ` .Emax` : upper energy limit (`::T`)
* ` .ΔE` : difference with respect to previous E used to calculate current E  (`::T`)
"""
function castInit(E::T, def::Def{T}) where T<:Real

    N = def.pos.N
    ℓ = def.orbit.ℓ

    nul = T(0)
    two = T(2)
    pot = def.potscr
  
    Nmin = getNmin(pot, 1:N)
    Emax = pot[N]
    Emin = pot[Nmin]
    if iszero(E)
        E = iszero(ℓ) ? 10Emax : 0.9Emin 
    else
        E = (Emin < E < Emax) ? E : iszero(ℓ) ? 10Emax : (Emin+Emax)/two
    end

    return Init(Emin, E, Emax, nul)
        
end

# ------------------------------------------------------------------------------
#                               init!{T}
# ------------------------------------------------------------------------------

function init!(init::Init{T}, ΔE::T, def::Def{T}) where T<:Real
    
    n′= def.orbit.n′     # radial quantum number (number of nodes)
    ℓ = def.orbit.ℓ
    nodes = def.pos.nodes
    
    two = T(2)

    Emin = init.Emin
    E = init.E
    Emax = init.Emax
        
    if ΔE ≠ T(0)
        if nodes > n′
            Emax = init.Emax = min(E, Emax)
            init.E = (E-ΔE+Emax)/two
        elseif nodes < n′
            Emin = init.Emin = max(E, Emin)
            init.E = (E-ΔE+Emin)/two
            println("Warning: code coverage: untested code path - analyze cause")
        else
            init.E = iseven(n′) ? E-ΔE : E+ΔE  # teken nu goed?
            init.ΔE = ΔE
        end
    else # node search:
        if nodes > n′
            Emax = init.Emax = min(E, Emax)
            init.E = ℓ > 0 ? (Emin+Emax)/two : init.E * T(1.1) 
        elseif nodes < n′
            Emin = init.Emin = max(E, Emin)
            init.E = (Emin+Emax)/two
        end
    end
        
    return init

end