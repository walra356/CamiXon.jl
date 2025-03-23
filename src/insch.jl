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
#                               insch.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#                       INSCH!(Z, E, grid, def)
# ------------------------------------------------------------------------------

@doc raw"""
    INSCH!(Z::Vector{Complex{T}}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real
 
Ansatz solution for the *inward* integration of the radial wave equation for the first ``k`` points 
on the [`CamiDiff.Grid`](@extref), where ``k`` is the Adams-Moulton order. The Ansatz is based on the WKB solution 
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
    k = grid.k
    Nuctp = def.pos.Nuctp
    pot = def.potscr

    p = Array{T,1}(undef,N)
    I = Array{T,1}(undef,N)
    P = Array{T,1}(undef,N)
    Q = Array{T,1}(undef,2k+1)
    
    for n=Nuctp:N
        p[n] = sqrt(abs(pot[n]-E))                         # quasi-classical momentum   
        I[n] = CamiDiff.grid_integration(p, grid, Nuctp:n)
        P[n] = exp(-I[n])/sqrt(p[n])                       # WKB solution
    end

    for n=Nuctp:N
        def.pos.Nb = P[n] > 1.0e-30 ? n : break
    end

    Nb = def.pos.Nb = def.pos.Nb - k
      
    Q = CamiDiff.grid_differentiation(P, grid, Nb-k:Nb+k)   # avoid lower end point correction by doubling range

    for n=Nb-k:Nb+k
        Z[n] = P[n] + im * Q[n-Nb+k+1]    
    end 
    
    return Z

end

# ------------------------------------------------------------------------------
#                       INSCH_WJ!(Z, E, grid, def) 
# kept for the record in "assets" folder
# ------------------------------------------------------------------------------
