# SPDX-License-Identifier: MIT

# author: Jook Walraven - 5-11-2024

# ==============================================================================
#                               pos.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#             Pos(Na, Nlctp, Nmin, Nuctp, Nb, N, nodes, ΔNuctp) 
# ------------------------------------------------------------------------------

@doc raw"""
    mutable struct Pos{T} where T<:Real

Type with fields:
* `   .Na::Int`: grid index of last leading point
* `.Nlctp::Int`: grid index of lower classical turning point
* ` .Nmin::Int`: grid index of (screened) potential minimum
* `.Nuctp::Int`: grid index of upper classical turning point
* `   .Nb::Int`: grid index first trailing point
* `    .N::Int`: grid index last point
* `.nodes::Int`: number of nodes in reduced wavefunction (r ≠ 0)
* `.ΔNlctp::Float64`: lctp offset with respect to Nlctp (1.0 ≤ ΔN ≤ 1.0)
* `.ΔNuctp::Float64`: uctp offset with respect to Nuctp (-1.0 ≤ ΔN ≤ 0.0)

Mutable struct to hold special grid indices as well as the number of nodes
and the (negative) offset of the exact uctp with respect to Nuctp.
`Pos` is one of the fields of the [`Def`](@ref) object
#### Examples:
```
julia> pos = Pos(1, 2, 3, 4, 5, 6, 7, 8.0, 9.0, 10.0);
julia> pos.Nuctp
4

julia> pos.Nuctp = 8;
julia> pos
Pos(1, 2, 3, 8, 5, 6, 7, 8.0, 9.0, 10.0)
```
"""
mutable struct Pos
    Na::Int
    Nlctp::Int
    Nmin::Int
    Nuctp::Int
    Nb::Int
    N::Int
    nodes::Int
    ΔNlctp::Float64
    ΔNuctp::Float64
    cWKB::Float64
end

# ------------------------------------------------------------------------------
#             castPos(E, Veff, grid) 
# ------------------------------------------------------------------------------

@doc raw"""
    castPos(E::T, Veff::Vector{T}, grid::Grid{T}) where T<:Real

Create the [`Pos`](@ref) object starting from the energy `E`, and the effective 
potential energy (screened Coulomb potential) `Veff[n]` tabulated on the [`Grid`](@ref). 
"""
function castPos(E::T, Veff::Vector{T}, grid::Grid{T}) where T<:Real

    N = grid.N
    k = grid.k

    firstabove = Veff[1] > E ? true : false
     lastabove = Veff[N] > E ? true : false

    lastabove || error("Error: Nuctp not found (resolve this by increasing Rmax)")

    Nmin = firstabove & lastabove ? getNmin(Veff, 1, N) : 1  
    # has minimum if (firstabove & lastabove == true)

    if firstabove & lastabove
        Nlctp = getNcut(E, Veff, 1, Nmin)
        Nuctp = getNcut(E, Veff, Nmin, N)
    else 
        Nlctp = 0
        Nuctp = getNcut(E, Veff, 1, N)
    end

    nodes = 0
    ΔNlctp = 0.0
    ΔNuctp = 0.0
    cWKB = 1e-7

    return Pos(k+1, Nlctp, Nmin, Nuctp, N-k, N, nodes, ΔNlctp, ΔNuctp, cWKB)
    
end

# ------------------------------------------------------------------------------
#             updatePos!(pos, E, Veff, grid) 
# ------------------------------------------------------------------------------

@doc raw"""
    updatePos!(pos::Pos, E::T, Veff::Vector{T}, grid::Grid{T}) where T<:Real

Update the [`Pos`](@ref) object starting from the energy `E`, and the effective 
potential energy (screened Coulomb potential) `Veff[n]` tabulated on the [`Grid`](@ref). 
"""
function updatePos!(pos::Pos, E::T, Veff::Vector{T}, grid::Grid{T}) where T<:Real

    N = grid.N
    k = grid.k

    firstabove = Veff[1] > E ? true : false
     lastabove = Veff[N] > E ? true : false

    lastabove || error("Error: Nuctp outside range (increase Rmax)")

    Nmin = (firstabove & lastabove) ? getNmin(Veff, 1, N) : 1  
    # Veff has minimum if (firstabove & lastabove) == true

    if firstabove & lastabove
        Nlctp = getNcut(E, Veff, 1, Nmin)
        Nuctp = getNcut(E, Veff, Nmin, N)
    else 
        Nlctp = 0
        Nuctp = getNcut(E, Veff, 1, N)
    end

    pos.Nlctp = Nlctp
    pos.Nmin = Nmin
    pos.Nuctp = Nuctp
    #pos.ΔNuctp = getΔNuctp(E, Veff, pos) 

    return pos
    
end

# ------------------------------------------------------------------------------
#                      getNmin(f, start, stop) 
# ------------------------------------------------------------------------------

@doc raw"""
    getNmin(f::Vector{T}, start::Int, stop:Int) where T<:Real
    getNmin(f::Vector{T}, itr::UnitRange) where T<:Real

Index corresponding to the absolute minimum of the discrete function ``f[n]`` truncated 
at the boundary of the interval ``start ≤ n ≤ stop``. 
Condition: ``f'[n]`` must be monotonically increasing or decreasing on the interval ``{start,stop}``.

NB. For a regular parabola the algorithm finds the index of the minimum. A truncated inverted parabola has 
two minima (at the boundaries of the interval). In this case the algorithm finds the index 
of the lowest of the two. If undecided the result is start.
"""
function getNmin(f::Vector{T}, start::Int, stop::Int) where T<:Real

     n = start
     m = stop 
    Δn = (m-n)÷2

    while Δn > 1
        if f[n+Δn] > f[n+Δn+1] # minimum to the right of n
            n += Δn 
        else                   # minimum to the left of n
            m -= Δn
        end
        Δn = (m-n)÷2
    end
    while f[n] > f[n+1]
        n += 1
        n < stop || break
    end
   
    return n
    
end
function getNmin(f::Vector{T}, itr::UnitRange) where T<:Real
    return getNmin(f, itr.start, itr.stop)
end

# ------------------------------------------------------------------------------
#                      getNmax(f, start, stop) 
# ------------------------------------------------------------------------------

@doc raw"""
    getNmax(f::Vector{T}, start::Int, stop:Int) where T<:Real
    getNmax(f::Vector{T}, itr::UnitRange) where T<:Real

Index corresponding to the absolute maximum of the discrete function ``f[n]`` truncated 
at the boundaries of the interval ``start ≤ n ≤ stop``. 
Condition: ``f'[n]`` must be monotonically increasing or decreasing on the interval ``{start,stop}``.

NB. For an inverted parabola the algorithm finds the index of the extremum. A regular parabola 
has two maxima (at the boundaries of the search interval). In this case the algorithm finds the 
index of the highest of the two. If undecided the result is start.
"""
function getNmax(f::Vector{T}, start::Int, stop::Int) where T<:Real

     n = start
     m = stop 
    Δn = (m-n)÷2

    while Δn > 1
        if f[n+Δn] < f[n+Δn+1] # maximum to the right of n
            n += Δn 
        else                   # maximum to the left of n
            m -= Δn
        end
        Δn = (m-n)÷2
    end
    while f[n] < f[n+1]
        n += 1
        n < stop || break
    end
   
    return n
    
end
function getNmax(f::Vector{T}, itr::UnitRange) where T<:Real
    return getNmax(f, itr.start, itr.stop)
end

# ------------------------------------------------------------------------------
#                      getNcut(val, f, start, stop) 
# ------------------------------------------------------------------------------

@doc raw"""
    getNcut(f0::T, f::Vector{T}, start::Int, stop::Int) where T<:Real

Index corresponding to the intersection point of the discrete function ``f[n]`` with the value ``f0`` 
in the interval ``start ≤ n ≤ stop``, 

```math
    f[n_cut] = f0.
```

Condition: ``f[n]`` must be monotonically increasing or decreasing on the interval ``start ≤ n ≤ stop``.

NB. For a monotonically *decreasing* function ``n_cut`` is approximated by the *largest* ``n`` for which ``f[n] > f0``.
For a monotonically *increasing* function  ``n_cut`` is approximated by the *smallest* ``n`` for which ``f[n] > f0``.
"""
function getNcut(f0::T, f::Vector{T}, start::Int, stop::Int) where T<:Real
    
     n = start
     m = stop    
    Δn = (m-n)÷2

    raising = (f[n] ≤ f0 ≤ f[m]) ? true : (f[n] ≥ f0 ≥ f[m]) ? false : error("Error: value outside function range")

    while Δn ≥ 1
        if raising
            if f[n+Δn] ≤ f0
                n += Δn   # raising ? n += Δn : n -= Δn
            elseif f[m-Δn] ≥ f0
                m -= Δn   # raising ? m -= Δn : m += Δn
            else 
                error("Error: search for intersection point undecided")
            end
        else
            if f[n+Δn] ≥ f0
                n += Δn   # raising ? n += Δn : n -= Δn
            elseif f[m-Δn] ≤ f0
                m -= Δn   # raising ? m -= Δn : m += Δn
            else 
                error("Error: search for intersection point undecided")
            end
        end
        Δn = (m-n)÷2
    end

    o = raising ? n+1 : n    # o is the index closest to val for which f[o] > val 

    return o

end
function getNcut(f0::T, f::Vector{T}, itr::UnitRange) where T<:Real
    getNcut(f0, f, itr.start, itr.stop)
end

# ------------------------------------------------------------------------------
#                      getΔNcut(f0, f, sense=fwd; ϵ = 1e-8, k = 7)
# ------------------------------------------------------------------------------

@doc raw"""
    getΔNcut(f0::T, f::Vector{T}, Ncut::Int, sense=fwd; ϵ = 1e-8, k = 7)

Offset of the exact intersection with respect to the index `Ncut` given as a Real number.
`ϵ` - convergence goal
`k` - order of a k+1 point Lagrange interpolation procedure based on a *linear* grid. 
Forward sense (`fwd`): value in the interval {0.0, 1.0}
Backward sense (`bwd`): value in the interval {-1.0, 0.0}
"""
function getΔNcut(f0::T, f::Vector{T}, Ncut::Int, sense=fwd; ϵ = 1e-8, k = 7) where T<:Real

    if CamiMath.isforward(sense)
        coords = CamiMath.lagrange_polynom(f, Ncut, Ncut+k, fwd)
        imin = 0.0   
        imax = 1.0   
        f1 = CamiMath.polynomial(coords, imax)
        f2 = CamiMath.polynomial(coords, imin) 
        f1 ≤ f0 ≤ f2 || error("Error: intersection condition 'f[Ncut] ≤ f0 ≤ f[Ncut+1]' violated")
    
        o = (imax+imin)/2.0   # forward offset w.r.t. Ncut
        while imax-imin > ϵ
            if CamiMath.polynomial(coords, o) ≤ f0
                imax = o
            elseif CamiMath.polynomial(coords, o) ≥ f0
                imin = o
            else
                println("noise level")
            end
            o = (imax+imin)/2.0
        end
    else    
        coords = CamiMath.lagrange_polynom(f, Ncut-k, Ncut, bwd)
        imin = -1.0
        imax = 0.0
        f1 = CamiMath.polynomial(coords, imin)
        f2 = CamiMath.polynomial(coords, imax)  
        f1 ≤ f0 ≤ f2 || error("Error: intersection condition 'f[Ncut-1] ≤ f0 ≤ f[Ncut]' violated")
    
        o = (imax+imin)/2.0   # backward offset w.r.t. Ncut
        while imax-imin > ϵ
            if CamiMath.polynomial(coords, o) ≥ f0
                imax = o
            elseif CamiMath.polynomial(coords, o) ≤ f0
                imin = o
            else
                println("noise level")
            end
            o = (imax+imin)/2.0
        end
    end

    return o
    
end

# ------------------------------------------------------------------------------
#                      getΔNuctp(E, Veff, pos) 
# ------------------------------------------------------------------------------

@doc raw"""
    getΔNuctp(E::T, Veff::Vector{T}, pos::Pos{T}) where T<:Real

Update the [`Pos`](@ref) object starting from the energy `E`, and the effective 
potential energy (screened Coulomb potential) `Veff[n]` tabulated on the [`Grid`](@ref). 
"""
function getΔNuctp(E::T, Veff::Vector{T}, pos::Pos) where T<:Real
    
    coords = CamiMath.lagrange_polynom(Veff, pos.Nuctp-3, pos.Nuctp, bwd)

    imax = 0.0
    imin = -1.0
    Δi = (imax+imin)/2.0   
    ϵ = 0.0001 

    while Δi  > ϵ
        if CamiMath.polynomial(coords, Δi) > E
            imax = Δi
        elseif CamiMath.polynomial(coords, Δi) < E
            imin = Δi
        end
        Δi = (imax+imin)/2.0
    end

    return Δi
    
end
