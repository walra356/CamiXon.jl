# ================== Pos(Na, Nlctp, Nmin, Nuctp, Nb, N, nodes) =================
@doc raw"""
    Pos(Na::Int, Nlctp::Int, Nmin::Int, Nuctp::Int, Nb::Int, N::Int, nodes::Int)

Type with fields:
* `   .Na`: grid index of last leading point (`::Int`)
* `.Nlctp`: grid index of classical turning point (`::Int`)
* ` .Nmin`: grid index of (screened) potential minimum (`::Int`)
* `.Nuctp`: grid index of classical turning point (`::Int`)
* `   .Nb`: grid index first trailing point (`::Int`)
* `    .N`: grid index last point (`::Int`)
* `.nodes`: number of nodes  (`::Int`)

Mutable struct to hold special grid indices as well as the number of nodes;
`Pos` is one of the fields of the [`Def`](@ref) object
#### Examples:
```
pos = Pos(1, 2, 3, 4, 5, 6, 7)
pos.Nuctp
 4

pos.Nuctp = 8
pos
 Pos(1, 2, 3, 8, 5, 6, 7)
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
end

# ===========   Grid (ID, name, Type, N, r, r′, h, r0, epn, epw, k) ============

@doc raw"""
    Def(T, atom, orbit, pot, scr, o1, o2, o3, pos, epn, k, am, matLD)

Type with fields:
* `    .T`: gridType (`::Type`)
* ` .atom`: atom object (`::Atom`)
* `.orbit`: orbit object (`::Orbit`)
* `  .pot`: tabulated potential function (`::Vector{T}`)
* `  .scr`: tabulated screening function (`::Vector{T}`)
* `   .o1`: vector of zero-filled matrices (`::Vector{Matrix{T}}`)
* `   .o2`: vector of zero-filled matrices (`::Vector{Matrix{T}}`)
* `   .o3`: vector of unit-filled matrices (`::Vector{Matrix{T}}`)
* `  .pos`: object containing Na, Nlctp, Nmin, Nuctp, Nb, N and nodes (`::Pos`)
* `  .epn`: number of endpoints trapezoidal correction - must be odd (`::Int`)
* `    .k`: Adams-Moulton order (`::Int`)
* `   .am`: Adams-Moulton weight coefficients (`::Vector{T}`)
* `.matLD`: Lagrangian differentiation matrix (`::Matrix{T}`)

The object `Def` is best created by the function [`castDef`](@ref).
"""
struct Def{T}
    T::Type
    atom::Atom
    orbit::Orbit
    pot::Vector{T}          # tabulated potential function
    scr::Vector{T}          # tabulated screening function
    o1::Vector{Matrix{T}}   # vector of zero-filled matrices
    o2::Vector{Matrix{T}}   # vector of zero-filled matrices
    o3::Vector{Matrix{T}}   # vector of unit-filled matrices
    pos::Pos                # object containing Nmin, Na, Nuctp, Nb, N and nodes
    epn::Int                # number of endpoints trapezoidal correction
    k::Int                  # Adams-Moulton order
    am::Vector{T}           # Adams-Moulton weight coefficients
    matLD::Matrix{T}        # Lagrangian differentiation matrix
end

# ================== castDef(grid, atom::Atom, orbit::Orbit) =================

function _defspecs(grid, atom, orbit)

    g = grid.name
    T = grid.T
    a = atom.name
    o = orbit.name

    return "Def created for $a $o on $g in $T"

end
# ..............................................................................
@doc raw"""
    castDef(grid::Grid{T}, atom::Atom, orbit::Orbit) where T <: Real

Create the [`Def`](@ref) object starting from the [`Grid`](@ref) object and the
atomic properties of the objects [`Atom`](@ref) and [`Orbit`](@ref).
#### Example:
```
atom = castAtom(Z=1, Q=0, M=1.00782503223, I=1//2, gI=5.585694713)
orbit = castOrbit(n=7, ℓ=2)
codata = castCodata(2018)
grid = autoGrid(atom, orbit, codata, Float64)
def = castDef(grid, atom, orbit);
    Atom created: Hydrogen - ¹H (Z = 1, Zc = 1, Q = 0, M = 1.00782503223, I = 1//2, gI = 5.585694713)
    Orbit created: 7d - (n = 7, n′ = 4, ℓ = 2)
    Grid created: exponential, Float64, Rmax = 207.0 (a.u.), Ntot = 400, h = 0.025, r0 = 0.00939821
    Def created for Hydrogen 7d on exponential grid in Float64
```
"""
function castDef(grid::Grid{T}, atom::Atom, orbit::Orbit; scr=nothing) where T <: Real
# ================================================================================
# castDef(grid, atom, orbit) # reference arrays
# ================================================================================
    N = grid.N
    r = grid.r
    k = grid.k
    epn = grid.epn
    Z = atom.Z
    ℓ = orbit.ℓ

    r[N]^(ℓ+1) < Inf || error("Error: numerical overflow (Inf)")

    Z = myconvert(T, Z)
    num = myconvert(T, ℓ*(ℓ + 1)//2)

    r1 = T(1.0e-100)  # quasi zero
    pot = ℓ > 0 ? [(-Z + num/r[n])/r[n] for n=1:N] : [-Z/r[n] for n=1:N]
    pot[1] = ℓ > 0 ? (-Z + num/r1)/r1 : -Z/r1
    pot = convert.(T,pot)
    scr = isnothing(scr) ? zeros(T,N) : scr
    o1 = [fill(myconvert(T,0), (2,2)) for n=1:N]
    o2 = [fill(myconvert(T,0), (2,2)) for n=1:N]
    o3 = [fill(myconvert(T,1), (2,2)) for n=1:N]
    pos = Pos(k+1, 0, 1, 0, N-k, N, 0)  # Pos(Na, Nlctp, Nmin, Nuctp, Nb, N, nodes)
    am = myconvert.(T, create_adams_moulton_weights(k; rationalize=true))
    matLD = myconvert.(T, create_lagrange_differentiation_matrix(k))

    println(_defspecs(grid, atom, orbit))

    return Def(T, atom, orbit, pot, scr, o1, o2, o3, pos, epn, k, am, matLD)

end

@doc raw"""
    initE(def::Def{T}; E=nothing) where T<:Real

Autogenerated seed value for the energy (default: no manual E seed)
#### Example:
```
codata = castCodata(2018)
atom = castAtom(Z=1, Q=0, M=1.00782503223, I=1//2, gI=5.585694713)
orbit = castOrbit(n=1, ℓ=0)
Ecal = convert(Float64, bohrformula(atom.Z, orbit.n))
grid = autoGrid(atom, orbit, codata; msg=false)
def = castDef(grid, atom, orbit)
  Atom created: Hydrogen - ¹H (Z = 1, Zc = 1, Q = 0, M = 1.00782503223, I = 1//2, gI = 5.585694713)
  Orbit created: 1s - (n = 1, n′ = 0, ℓ = 0)

E = initE(def); println("E = $E")
  E = -0.03508495857961283

E = initE(def; E=Ecal); println("E = $E")
  E = -0.5
```
"""
function initE(def::Def{T}; E=nothing) where T<:Real

    isnothing(E) || return myconvert(T, E)

    N = def.pos.N
    ℓ = def.orbit.ℓ
    v = def.pot
    s = def.scr

    pot = v .+ s

    Emax = pot[N]
    Emin = minimum(pot[2:N])

    E = iszero(ℓ) ? 2.0Emax : 0.9Emin

    return E

end

# ============================== get_Na(Z, def) ================================
@doc raw"""
    get_Na(Z::Vector{Complex{T}}, def::Def{T}) where T<:Real

Grid index of the starting point for outward numerical integration. This is
the first point where the integration threshold value (1.0e-10) is exceeded.
"""
function get_Na(Z::Vector{Complex{T}}, def::Def{T}) where T<:Real
# ==============================================================================
#  grid index of starting point outward numerical integration
# ==============================================================================
    k = def.k

    ref = T(1.0e-10)

    Na = findfirst(x -> abs(x) > ref, real(Z))
    Na = isnothing(Na) ? k+1 : Na > 0 ? max(k+1, Na) : k+1

    return Na

end

# ============================= get_Nb(Z, def) =================================
@doc raw"""
    get_Nb(Z::Vector{Complex{T}}, def::Def{T}) where T<:Real

Grid index of the stopping for outward numerical integration. This is
the last point where the integration threshold value (1.0e-10) is exceeded.
"""
function get_Nb(Z::Vector{Complex{T}}, def::Def{T}) where T<:Real
# ==============================================================================
#  grid index of starting point inward numerical integration
# ==============================================================================
    k = def.k
    N = def.pos.N

    ref = T(1.0e-10)

    Nb = findlast(x -> abs(x) > ref, real(Z))
    Nb = isnothing(Nb) ? N-k : Nb > 0 ? min(N-k, Nb) : N-k

    return Nb

end
# ============================= get_Nb(Z, def) =================================
@doc raw"""
    get_Nmin(def::Def{T}) where T<:Real

Grid index of the minimum of the screened potential curve. By definition
`get_Nmin(def) = 1` for zero orbital angular momentum (``ℓ=0``).
"""
function get_Nmin(def::Def{T}) where T<:Real
# ==============================================================================
#  grid index of potential minimum
# ==============================================================================
    N = def.pos.N
    v = def.pot
    s = def.scr

    pot = v .+ s

    n = 1
    while pot[n+1] < pot[n]
        n < N-1 ? n += 1 : break
    end

    return n

end

# ============================= get_Nuctp(E, def) ==============================
@doc raw"""
    get_Nlctp(E::T, def::Def{T}) where T<:Real

Grid index of the *lower classical turning point * of the screened potential
curve. By definition `get_Nlctp(E, def) = 2` for zero orbital angular
momentum (``ℓ=0``).
"""
function get_Nlctp(E::T, def::Def{T}) where T<:Real
# ==============================================================================
#  grid index of lower classical turning point
# ==============================================================================
    Nmin = def.pos.Nmin
    N = def.pos.N
    v = def.pot
    s = def.scr
    ℓ = def.orbit.ℓ

    ℓ > 0 || return 0

    pot = v .+ s

    E < pot[Nmin] && println("Warning: E < Emin - force Nuctp = Nmin")

    n = 2

    while E < pot[n]     # below classical threshhold
        n < Nmin ? n += 1 : break
    end

    n < N || error("Error: inner classical turning point outside grid")

    return n

end

# ============================= get_Nuctp(E, def) ==============================
@doc raw"""
    get_Nuctp(E::T, def::Def{T}) where T<:Real

Grid index of the *upper classical turning point* of the screened potential
curve. By definition `get_Nuctp(E, def) = N-1` for zero orbital angular
momentum (``ℓ=0``).
"""
function get_Nuctp(E::T, def::Def{T}) where T<:Real
# ==============================================================================
#  grid index of upper classical turning point
# ==============================================================================
    Nmin = def.pos.Nmin
    N = def.pos.N
    v = def.pot
    s = def.scr

    pot = v .+ s

    n = N-1

    E < pot[Nmin] && println("Warning: E < Emin - force Nuctp = Nmin")

    while E < pot[n]     # below classical threshhold
        n > Nmin ? n -= 1 : break
    end

    n < N || error("Error: (outer) classical turning point outside grid")

    return n

end

# ======================= get_nodes(Z, def) ==============================

@doc raw"""
    get_nodes(Z::Vector{Complex{T}}, def::Def{T}) where T<:Real

Number of nodes (excluding the origin) of the reduced radial wavefunction
χ(r) = real(Z).
"""
function get_nodes(Z::Vector{Complex{T}}, def::Def{T}) where T<:Real

    Na = def.pos.Na
    Nb = def.pos.Nb
    c = [sign(real(Z[n])*real(Z[n+2])) for n=Na:2:Nb-2]
    c = findall(x -> x == -1, c)

    nodes = length(c)

    return nodes

end
