# ================== Pos(Na, Nlctp, Nmin, Nuctp, Nb, N, nodes) =================
@doc raw"""
    Pos(Na::Int, Nlctp::Int, Nmin::Int, Nuctp::Int, Nb::Int, N::Int, nodes::Int, cWKB::Float64)

Type with fields:
* `   .Na`: grid index of last leading point (`::Int`)
* `.Nlctp`: grid index of lower classical turning point (`::Int`)
* ` .Nmin`: grid index of (screened) potential minimum (`::Int`)
* `.Nuctp`: grid index of upper classical turning point (`::Int`)
* `   .Nb`: grid index first trailing point (`::Int`)
* `    .N`: grid index last point (`::Int`)
* `.nodes`: number of nodes  (`::Int`)
* ` .cWKB`: WKB threshold level determining Na and Nb (`::Float64`)

Mutable struct to hold special grid indices as well as the number of nodes;
`Pos` is one of the fields of the [`Def`](@ref) object
#### Examples:
```
pos = Pos(1, 2, 3, 4, 5, 6, 7, 8)
pos.Nuctp
    4

pos.Nuctp = 8
pos
    Pos(1, 2, 3, 9, 5, 6, 7, 8)
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
    cWKB::Float64
end

# ===========   Grid (ID, name, Type, N, r, r′, h, r0, epn, epw, k) ============

@doc raw"""
    Def(T, atom, orbit, pot, scr, o1, o2, o3, pos, epn, k, am, matLD)

Type with fields:
* `     .T`: gridType (`::Type`)
* `  .atom`: atom object (`::Atom`)
* ` .orbit`: orbit object (`::Orbit`)
* `.codata`: codata object (`::Codata`)
* `   .pot`: tabulated potential function (`::Vector{T}`)
* `   .scr`: tabulated screening function (`::Vector{T}`)
* `    .o1`: vector of zero-filled matrices (`::Vector{Matrix{T}}`)
* `    .o2`: vector of zero-filled matrices (`::Vector{Matrix{T}}`)
* `    .o3`: vector of unit-filled matrices (`::Vector{Matrix{T}}`)
* `   .pos`: object containing Na, Nlctp, Nmin, Nuctp, Nb, N and nodes (`::Pos`)
* `   .epn`: number of endpoints trapezoidal correction - must be odd (`::Int`)
* `     .k`: Adams-Moulton order (`::Int`)
* `    .am`: Adams-Moulton weight coefficients (`::Vector{T}`)
* ` .matLD`: Lagrangian differentiation matrix (`::Matrix{T}`)

The object `Def` is best created with the function [`castDef`](@ref).
"""
struct Def{T}
    T::Type
    atom::Atom
    orbit::Orbit
    codata::Codata
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

# ========= castDef(grid, atom::Atom, orbit::Orbit, codata::Codata) ============

function _defspecs(grid, atom, orbit)

    g = grid.name
    T = grid.T
    a = atom.element.name
    o = orbit.name

    return "Def created for $a $o on $g grid of $(grid.N) points"

end
# ..............................................................................
@doc raw"""
    castDef(grid::Grid{T}, atom::Atom, orbit::Orbit, codata::Codata[; scr=nothing[, msg=true]]) where T <: Real

Create the [`Def`](@ref) object starting from the [`Grid`](@ref) object and the
atomic properties of the objects [`Atom`](@ref) and [`Orbit`](@ref).
Optional: scr (supply screening array)
#### Example:
```
codata = castCodata(2018)
atom = castAtom(Z=1, A=1, Q=0)
orbit = castOrbit(n=7, ℓ=2)
grid = autoGrid(atom, orbit, Float64)
def = castDef(grid, atom, orbit, codata);
    Element created: H, hydrogen, Z=1, weight=1.008
    Isotope created: ¹H, hydrogen, Z=1, A=1, N=0, R=0.8783, M=1.007825032, I=1/2⁺, μI=2.792847351, Q=0.0, RA=99.9855%, (stable)
    Atom created: hydrogen, neutral atom, ¹H, Z=1, A=1, Q=0, Zc=1
    Orbital: 7d
        principal quantum number: n = 7
        radial quantum number: n′ = 4 (number of nodes in radial wavefunction)
        orbital angular momentum of valence electron: ℓ = 2
    Grid created: exponential, Float64, Rmax = 207.0 a.u., Ntot = 400, h = 0.025, r0 = 0.00939821
    Def created for hydrogen 7d on exponential grid of 400 points
```
"""
function castDef(grid::Grid{T}, atom::Atom, orbit::Orbit, codata::Codata; scr=nothing, msg=true) where T <: Real
# ================================================================================
# castDef(grid, atom, orbit, codata) # reference arrays
# ================================================================================
    N = grid.N
    r = grid.r
    k = grid.k
    epn = grid.epn
    Z = atom.Z
    ℓ = orbit.ℓ

    r[N]^(ℓ+1) < Inf || error("Error: numerical overflow (r[N]^(ℓ+1) -> Inf)")

    Z = convert(T, Z)
    num = convert(T, ℓ*(ℓ + 1)//2)

    r1 = T(1.0e-100)  # quasi zero
    pot = ℓ > 0 ? [(-Z + num/r[n])/r[n] for n=1:N] : [-Z/r[n] for n=1:N]
    pot[1] = ℓ > 0 ? (-Z + num/r1)/r1 : -Z/r1
    pot = convert.(T,pot)
    scr = isnothing(scr) ? zeros(T,N) : scr
    o1 = [fill(convert(T,0), (2,2)) for n=1:N]
    o2 = [fill(convert(T,0), (2,2)) for n=1:N]
    o3 = [fill(convert(T,1), (2,2)) for n=1:N]
    pos = Pos(k+1, 0, 1, 0, N-k, N, 0, 1.0e-7)  # Pos(Na, Nlctp, Nmin, Nuctp, Nb, N, nodes)
    am = convert.(T, create_adams_moulton_weights(k; rationalize=true))
    matLD = convert.(T, create_lagrange_differentiation_matrix(k))

    msg && println(_defspecs(grid, atom, orbit))

    return Def(T, atom, orbit, codata, pot, scr, o1, o2, o3, pos, epn, k, am, matLD)

end

@doc raw"""
    initE(def::Def{T}) where T<:Real

Autogenerated seed value for the energy
#### Example:
```
codata = castCodata(2018)
atom = castAtom(Z=1, A=1, Q=0; msg=false)
orbit = castOrbit(n=1, ℓ=0; msg=false)
grid = autoGrid(atom, orbit, Float64; msg=false)
def = castDef(grid, atom, orbit, codata);
    Def created for hydrogen 1s on exponential grid of 100 points

E = initE(def); println("E = $E")
    E = -0.03508495857961283
```
"""
function initE(def::Def{T}) where T<:Real

    N = def.pos.N
    ℓ = def.orbit.ℓ
    v = def.pot
    s = def.scr

    pot = v .+ s

    Emax = pot[N]
    Emin = minimum(pot[2:N])

    E = iszero(ℓ) ? 10.0Emax : 0.9Emin

    return E

end

# ============================== get_Na(Z, def) ================================
@doc raw"""
    get_Na(Z::Vector{Complex{T}}, def::Def{T}) where T<:Real

Grid index of the starting point for outward numerical integration. This is
`k+1` or the point marking the end of the quasiclassical region below the
lower classical turning point (`lctp`) as marked by the WKB threshold value
(`def.pos.cWKB`).
#### Example:
```
Ecal, grid, def, adams = demo_hydrogen(n=1, ℓ=0)
E, def, adams, Z = adams_moulton_master(E, codata, grid, def, adams; Δν=Value(1,"kHz"), imax=25, msg=false);
    Orbital: 1s
        principal quantum number: n = 1
        radial quantum number: n′ = 0 (number of nodes in radial wavefunction)
        orbital angular momentum of valence electron: ℓ = 0
    Grid created: exponential, Float64, Rmax = 63.0 a.u., Ntot = 100, h = 0.1, r0 = 0.00286033
    Def created for hydrogen 1s on exponential grid

Na = get_Na(Z, def)
println("k + 1 = $(grid.k+1); Na = $Na")
    k + 1 = 8; Na = 8

Na == def.pos.Na
    true
```
"""
function get_Na(Z::Vector{Complex{T}}, def::Def{T}) where T<:Real
# ==============================================================================
#  grid index of starting point outward numerical integration
# ==============================================================================
    k = def.k
    cWKB = def.pos.cWKB
    #ref = T(1.0e-10)

    Na = findfirst(x -> abs(x) > cWKB, real(Z))
    Na = isnothing(Na) ? k+1 : Na > 0 ? max(k+1, Na) : k+1

    return Na

end

# ============================= get_Nb(Z, def) =================================
@doc raw"""
    get_Nb(Z::Vector{Complex{T}}, def::Def{T}) where T<:Real

Grid index of the stopping for outward numerical integration. This is `N-k-1`
or the point marking the start of the quasiclassical region above the
upper classical turning point (`Nuctp`) as marked by the WKB threshold value
(`def.pos.cWKB`).
#### Example:
```
Ecal, grid, def, adams = demo_hydrogen(n=1, ℓ=0)
E, def, adams, Z = adams_moulton_master(E, codata, grid, def, adams; Δν=Value(1,"kHz"), imax=25, msg=false);
    Orbital: 1s
        principal quantum number: n = 1
        radial quantum number: n′ = 0 (number of nodes in radial wavefunction)
        orbital angular momentum of valence electron: ℓ = 0
    Grid created: exponential, Float64, Rmax = 63.0 a.u., Ntot = 100, h = 0.1, r0 = 0.00286033
    Def created for hydrogen 1s on exponential grid

Nb = get_Nb(Z, def)
println("N - k - 1 = $(grid.N-grid.k-1); Nb = $Nb")
    N - k - 1 = 92; Nb = 92

Nb == def.pos.Nb
    true
```
"""
function get_Nb(Z::Vector{Complex{T}}, def::Def{T}) where T<:Real
# ==============================================================================
#  grid index of starting point inward numerical integration
# ==============================================================================
    k = def.k
    N = def.pos.N
    cWKB = def.pos.cWKB

    cWKB = T(cWKB)

    Nb = findlast(x -> abs(x) > cWKB, real(Z))
    Nb = isnothing(Nb) ? N-k : Nb > 0 ? min(N-k, Nb) : N-k

    return Nb

end
# ============================= get_Nb(Z, def) =================================
@doc raw"""
    get_Nmin(def::Def{T}) where T<:Real

Grid index of the minimum of the screened potential curve. By definition
`get_Nmin(def) = 1` for zero orbital angular momentum (`ℓ=0`).
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
momentum (`ℓ=0`).
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

# ======================= count_nodes(Z, def) ==============================

@doc raw"""
    count_nodes(Z::Vector{Complex{T}}, def::Def{T}) where T<:Real

Number of nodes (excluding the origin) of the reduced radial wavefunction
χ(r) = real(Z).

#### Example:
```
atom = castAtom(Z=1, A=1, Q=0, msg=false);
orbit = castOrbit(n=3, ℓ=2, msg=false);
grid = autoGrid(atom, orbit, Float64; Nboost=1, epn=5, k=7, msg=false);
def = castDef(grid.T, atom, orbit, codata);
    Def created for hydrogen 3d on exponential grid of 200 points

E = convert(setT, bohrformula(atom.Z, orbit.n));
adams = castAdams(E, grid, def);
E, def, adams, Z = adams_moulton_master(E, codata, grid, def, adams; Δν=Value(1,"kHz"), imax=25, msg=false);

o = count_nodes(Z, def); println("node count: $o nodes")
    node count: 0 nodes
```
"""
function count_nodes(Z::Vector{Complex{T}}, def::Def{T}) where T<:Real

    Na = def.pos.Na
    Nuctp = def.pos.Nuctp
    nodes = 0
    for n=Na+1:Nuctp-1
        if real(Z[n-1])*real(Z[n]) < 0
            nodes += 1
        end
    end

    return nodes

end
