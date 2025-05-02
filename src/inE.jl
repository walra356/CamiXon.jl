# SPDX-License-Identifier: MIT

# Copyright (c) 2025 Jook Walraven <69215586+walra356@users.noreply.github.com> and contributors

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# ==============================================================================
#                               inE.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#                               InE{T}
# ------------------------------------------------------------------------------

@doc raw"""
    InE{T} where T<:Real

Energy object with fields:

* ` .Emin` : lower energy limit (`::T`)
* ` .E` : trial energy (`::T`)
* ` .Emax` : upper energy limit (`::T`)
* ` .ΔE` : difference with respect to previous E used to calculate current E  (`::T`)
"""
mutable struct InE{T}

    Emin::T     # lowest allowed energy
    E::T        # trial energy
    Emax::T     # practical highest energy
    ΔE::T       # ΔE with respect to previous E to calculate E trial in use
    
end

# ------------------------------------------------------------------------------
#                               castInE{T}
# ------------------------------------------------------------------------------

@doc raw"""
    castInE(E::T, def::Def{T}) where T<:Real

Create and initialize the energy object [`InE`](@ref) according to the 
specifications given in the [`Def`](@ref) object.

* ` .Emin` : lower energy limit (`::T`)
* ` .E` : trial energy (`::T`)
* ` .Emax` : upper energy limit (`::T`)
* ` .ΔE` : difference with respect to previous E used to calculate current E  (`::T`)
"""
function castInE(E::T, def::Def{T}) where T<:Real

    N = def.pos.N
    ℓ = def.spinorbit.ℓ

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

    return InE(Emin, E, Emax, nul)
        
end

# ------------------------------------------------------------------------------
#                               inE!{T}
# ------------------------------------------------------------------------------

function inE!(inE::InE{T}, ΔE::T, def::Def{T}) where T<:Real
    
    n′= def.spinorbit.n′     # radial quantum number (number of nodes)
    ℓ = def.spinorbit.ℓ
    nodes = def.pos.nodes
    
    two = T(2)

    Emin = inE.Emin
    E = inE.E
    Emax = inE.Emax
        
    if ΔE ≠ T(0)
        if nodes > n′
            Emax = inE.Emax = min(E, Emax)
            inE.E = (E-ΔE+Emax)/two
        elseif nodes < n′
            Emin = inE.Emin = max(E, Emin)
            inE.E = (E-ΔE+Emin)/two
            println("Warning: code coverage: untested code path - analyze cause")
        else
            inE.E = iseven(n′) ? E-ΔE : E+ΔE  # teken nu goed?
            inE.ΔE = ΔE
        end
    else # node search:
        if nodes > n′
            Emax = inE.Emax = min(E, Emax)
            inE.E = ℓ > 0 ? (Emin+Emax)/two : inE.E * T(1.1) 
        elseif nodes < n′
            Emin = inE.Emin = max(E, Emin)
            inE.E = (Emin+Emax)/two
        end
    end
        
    return inE

end