# Grid

The [`CamiDiff.Grid`](@extref) object is the backbone for the numerical procedure on a non-uniform
grid. Its principal fields are `grid.r` and `grid.r′`, which are discrete
functions of `N` elements representing the grid function and its derivative.

#### autoGrid
```@docs
autoGrid(atom::Atom, orbit::Orbit, T::Type; p=0, rmax=0, N=0, polynom=[], epn=5, k=5, msg=false)
```

#### autoNtot
```@docs
autoNtot(n::Int)
```
#### autoRmax
```@docs
autoRmax(atom::Atom, n::Int; rmax=0.0)
```
#### autoPrecision
```@docs
autoPrecision(rmax::T, ℓ = 0) where T<:Real
```

## Def

The `Def` object serves to define the problem to be solved and to contain in
the field `def.Z` the solution as a discrete function of `N` elements.

```@docs
Def{T}
```
#### castDef
```@docs
castDef(grid::CamiDiff.Grid{T}, atom::Atom, spinorbit::Spinorbit, codata::Codata; pos=nothing, scr=nothing, msg=false) where T <: Real
```

## Pos

The `Pos` object serves within [`Def`](@ref) object to contain the position
indices `def.Na`, `def.Nb`, `def.Nlctp`, `def.Nmin`, `def.Nuctp` used in
Adams-Moulton integration. These positions are contained in the fields
`def.pos.Na`, `def.pos.Nb`, `def.pos.Nlctp`, `def.pos.Nmin`, `def.pos.Nuctp`.

```@docs
Pos
```
#### castPos
```@docs
castPos(E::T, Veff::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
```
#### updatePos!
```@docs
updatePos!(pos::Pos, E::T, Veff::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
```
#### listPos
```@docs
listPos(pos::Pos; msg=true)
```
#### Pos-related functions
```@docs
getNmin(f::Vector{T}, start::Int, stop::Int) where T<:Real
getNmax(f::Vector{T}, start::Int, stop::Int) where T<:Real
getNcut(f0::T, f::Vector{T}, start::Int, stop::Int) where T<:Real
getΔNcut(f0::T, f::Vector{T}, Ncut::Int, sense=fwd; ϵ = 1e-8, k = 7) where T<:Real
```