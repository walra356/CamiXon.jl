# ======================== Pos(Nmin, Na, Nctp, Nb, N, nodes) ===================
"""
    Pos(Nmin::Int, Na::Int, Nctp::Int, Nb::Int, N::Int, nodes::Int)

Type with fields:
* ` .Nmin`:grid index of (screened) potential minimum
* `   .Na`:grid index of last leading point
* ` .Nctp`:grid index of classical turning point
* `   .Nb`:grid index first trailing point
* `    .N`:grid index last point
* `.nodes`:number of nodes

Mutable struct to hold special grid indices as well as the number of nodes
### Example:
```
pos = Pos(1,2,3,4,5,6)
pos.Nctp
 3

pos.Nctp = 7
pos
 Pos(1, 2, 7, 4, 5, 6)
```
"""
mutable struct Pos
    Nmin::Int
    Na::Int
    Nctp::Int
    Nb::Int
    N::Int
    nodes::Int
end

# ===========   Grid (ID, name, Type, N, r, r′, h, r0, epn, epw, k) ============

"""
    Def(T, atom, orbit, pot, scr, o1, o2, o3, pos, k, am, matLD)

Type with fields:
* `    .T`::Type--------------- gridType
* ` .atom`::Atom                atom object
* `.orbit`::Orbit               orbit object
* `  .pot`::Vector{T}           tabulated potential function
* `  .scr`::Vector{T}           tabulated screening function
* `   .o1`::Vector{Matrix{T}}   vector of zero-filled matrices
* `   .o2`::Vector{Matrix{T}}   vector of zero-filled matrices
* `   .o3`::Vector{Matrix{T}}   vector of unit-filled matrices
* `  .pos`::Pos                 pos object for Nmin, Na, Nctp, Nb, N and nodes
* `    .k`::Int                 Adams-Moulton order
* `   .am`::Vector{T}           Adams-Moulton weight coefficients
* `.matLD`::Matrix{T}           Lagrangian differentiation matrix

The type `Grid` is best created by the function `createGrid`.
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
    pos::Vector{Int}        # position indices of 4 grid points [N, Na, Nb, Nctp]
    k::Int                  # Adams-Moulton order
    am::Vector{T}           # Adams-Moulton weight coefficients
    matLD::Matrix{T}        # Lagrangian differentiation matrix
end

"""
    createDef(grid::Grid{T}, atom::Atom, orbit::Orbit) where T <: Real

Create the Def object starting from the Grid and atomic properties.
"""
function createDef(grid::Grid{T}, atom::Atom, orbit::Orbit) where T <: Real
# ================================================================================
# createDef(grid, atom, orbit) # reference arrays
# ================================================================================
    N = grid.N
    r = grid.r
    k = grid.k
    Z = atom.Z
    ℓ = orbit.ℓ

    r[N]^(ℓ+1) < Inf || error("Error: numeric overflow, consider arbitrary precision type)")

    Z = myconvert(T, Z)
    num = myconvert(T, ℓ*(ℓ + 1)//2)

    r1 = T(1.0e-100)  # quasi zero
    pot = ℓ > 0 ? [(-Z + num/r[n])/r[n] for n=1:N] : [-Z/r[n] for n=1:N]
    pot[1] = ℓ > 0 ? (-Z + num/r1)/r1 : -Z/r1
    pot = convert.(T,pot)
    scr = zeros(T,N)
    o1 = [fill(myconvert(T,0), (2,2)) for n=1:N]
    o2 = [fill(myconvert(T,0), (2,2)) for n=1:N]
    o3 = [fill(myconvert(T,1), (2,2)) for n=1:N]
    pos = Pos(1, 1, 0, N, N, 0)  # placeholders for Nmin, Na, Nctp, Nb, N and nodes
    am = myconvert.(T, create_adams_moulton_weights(k; rationalize=true))
    matLD = myconvert.(T, create_lagrange_differentiation_matrix(k))

    return Def(T, atom, orbit, pot, scr, o1, o2, o3, pos, k, am, matLD)

end
