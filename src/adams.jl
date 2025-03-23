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

end

# ------------------------------------------------------------------------------
#                   castAdams(E, grid, def)
# ------------------------------------------------------------------------------

@doc raw"""
    castAdams(E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real

Initiates the [`Adams`](@ref) object.
"""
function castAdams(E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real

    castPos(E, def.potscr, grid)

    G = matG(E, grid, def)
    σ = matσ(E, grid, def)
    M = matMinv(E, grid, def)

    return Adams(G, σ, M)

end

# ------------------------------------------------------------------------------
#                   updateAdams!(adams, E, grid, def)
# ------------------------------------------------------------------------------

@doc raw"""
    updateAdams!(adams::Adams{T}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real

Update [`Adams`](@ref) object.
"""
function updateAdams!(adams::Adams{T}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real

G = matG(E, grid, def)
σ = matσ(E, grid, def)
M = matMinv(E, grid, def)

return Adams(G, σ, M)

end
