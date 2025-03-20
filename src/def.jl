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
#                               def.jl
# ==============================================================================

# =========== Def(atom, orbit, codata, pot, scr, potscr, G, σ, Minv, pos, epn, k, am, matLD) ============

@doc raw"""
    Def(T, atom, spinorbit, pot, scr, o1, o2, o3, pos, epn, k, am, matLD)

Type with fields:
* `  .atom::Atom`             : atom object
* ` .spinorbit::Spinorbit`    : spinorbit object
* `.codata::Codata`           : codata object
* `   .pot::Vector{T}`        : tabulated potential function
* `   .scr::Vector{T}`        : tabulated screening function
* `.potscr::Vector{T}`        : tabulated screened potential function
* `     .G::Vector{Matrix{T}}`: vector of zero-filled matrices
* `     .σ::Vector{Matrix{T}}`: vector of zero-filled matrices
* `  .Minv::Vector{Matrix{T}}`: vector of zero-filled matrices
* `   .pos::Pos`              : object with fields Na, Nlctp, Nmin, Nuctp, Nb, N and nodes
* `   .epn::Int`              : number of endpoints trapezoidal correction - must be odd
* `     .k::Int`              : Adams-Moulton order 
* `    .am::Vector{T}`        : Adams-Moulton weight coefficients
* ` .matLD::Matrix{T}`        : Lagrangian differentiation matrix

The object `Def` is best created with the function [`castDef`](@ref).
"""
struct Def{T}
    atom::Atom
    spinorbit::Spinorbit
    codata::Codata
    pot::Vector{T}          # tabulated potential function
    scr::Vector{T}          # tabulated screening function
    potscr::Vector{T}       # tabulated screened potential function
    G::Vector{Matrix{T}}    # uninitialized vector of matrices
    σ::Vector{Matrix{T}}    # uninitialized vector of matrices
    Minv::Vector{Matrix{T}} # uninitialized vector of matrices
    pos::Pos                # object containing Nmin, Na, Nuctp, Nb, N and nodes
    epn::Int                # number of endpoints trapezoidal correction
    k::Int                  # Adams-Moulton order
    am::Vector{T}           # Adams-Moulton weight coefficients
    matLD::Matrix{T}        # Lagrangian differentiation matrix
end




# ========= castDef(grid, atom::Atom, orbit::Orbit, codata::Codata) ============

function _defspecs(grid, atom, spinorbit)

    g = grid.name
    o = spinorbit.name
    Q = atom.Q

    strQ = abs(Q) > 1 ? sup(abs(Q)) : ""
    strQ = Q > 0 ? (strQ * 'ᐩ') : Q < 0 ? (strQ * 'ᐨ') : ""
    strN = Q ≠ 0 ? " ion" : ", neutral atom"

    str = atom.isotope.symbol * strQ

    return "Def created for "* str *":$o on $g grid of $(grid.N) points"

end
function _defspecs(grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real

    g = grid.name
    o = def.spinorbit.name
    Q = def.atom.Q

    strQ = abs(Q) > 1 ? sup(abs(Q)) : ""
    strQ = Q > 0 ? (strQ * 'ᐩ') : Q < 0 ? (strQ * 'ᐨ') : ""
    strN = Q ≠ 0 ? " ion" : ", neutral atom"

    str = def.atom.isotope.symbol * strQ

    return str * ":$o on $g grid of $(grid.N) points"

end
# ..............................................................................
@doc raw"""
    castDef(grid::CamiDiff.Grid{T}, atom::Atom, spinorbit::Spinorbit, codata::Codata [; pos=nothing, [scr=nothing[, msg=false]]) where T <: Real

Create the [`Def`](@ref) object starting from the [`CamiDiff.Grid`](@extref) object and the
atomic properties of the objects [`Atom`](@ref) and [`Orbit`](@ref).
Optional: scr (supply screening array)
#### Example:
```
julia> codata = castCodata(2018)
julia> atom = castAtom(Z=1, A=1, Q=0);
julia> orbit = castOrbit(n=7, ℓ=2);
julia> grid = autoGrid(atom, orbit, Float64);

julia> castDef(grid, atom, orbit, codata, msg=true);
Def created for hydrogen 7d on exponential grid of 400 points
```
"""
function castDef(grid::CamiDiff.Grid{T}, atom::Atom, spinorbit::Spinorbit, codata::Codata; pos=nothing, scr=nothing, msg=false) where T <: Real
    # ================================================================================
    # castDef(grid, atom, orbit, codata) # reference arrays
    # ================================================================================
        N = grid.N
        r = grid.r
        k = grid.k
        epn = grid.epn
        ℓ = spinorbit.ℓ
    
        #r[N]^(ℓ+1) < Inf || error("Error: numerical overflow (r[N]^(ℓ+1) -> Inf)")
    
        Z = convert(T, atom.Z)
        num = convert(T, ℓ*(ℓ + 1)//2)
    
        r1 = T(eps(Float64))  # quasi zero
        pot = ℓ > 0 ? [(-Z + num/r[n]) for n=1:N]./r : [-Z/r[n] for n=1:N]
        pot[1] = ℓ > 0 ? (-Z + num/r1)/r1 : -Z/r1
        pot = convert.(T,pot)
        scr = isnothing(scr) ? zeros(T,N) : convert.(T,scr)
        potscr = pot .+ scr
        G = [fill(convert(T,0), (2,2)) for n=1:N]
        σ = [fill(convert(T,0), (2,2)) for n=1:N]
        Minv = [fill(convert(T,1), (2,2)) for n=1:N]
        pos = isnothing(pos) ? Pos(k+1, 0, 1, 0, N-k, N, 0, 0.0, 0.0, 1e-7)  : pos # Pos(Na, Nlctp, Nmin, Nuctp, Nb, N, nodes, ΔNlctp, ΔNuctp, cWKB)
        am = convert.(T, CamiDiff.create_adams_moulton_weights(k; rationalize=true))
        matLD = convert.(T, CamiDiff.create_lagrange_differentiation_matrix(k))
    
        msg && println(_defspecs(grid, atom, spinorbit))
    
        return Def(atom, spinorbit, codata, pot, scr, potscr, G, σ, Minv, pos, epn, k, am, matLD)
    
    end

