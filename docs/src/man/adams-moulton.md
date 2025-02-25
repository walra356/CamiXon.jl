# Adams-Moulton integration

The Adams-Moulton method is used for numerical integration of the reduces
radial wave equation. In the present implementation it is constructed on top
the objects [`Atom`](@ref), [`Orbit`](@ref), [`CamiDiff.Grid`](@extref), [`Def`](@ref)
and [`Adams`](@ref) using 5 globally defined instances called `atom`, `orbit`,
`grid`, `def` and `adams`.

## Adams

The `Adams` object serves to hold the Adams-Moulton integration matrices
`matG`, `matσ`, `matMinv` as well as the *actual* normalized solution `Z` in
the form of a tabulated function of `N` elements.

```@docs
Adams
castAdams(E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real
updateAdams!(adams::Adams{T}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real
```

#### Adams related functions

```@docs
matG(E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real
matσ(E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real
matMinv(E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real
```

## Adams-Moulton numerical solution of the radial wave equation
```@docs
Init{T} where T<:Real
adams_moulton_solve!(Z::Vector{Complex{T}}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real
adams_moulton_solve_refine!(Z::Vector{Complex{T}}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real
```

## Radial integration - outward
```@docs
OUTSCH!(Z::Vector{Complex{T}}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real
OUTSCH_WJ!(Z::Vector{Complex{T}}, grid::CamiDiff.Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real
OUTSCH_WKB!(Z::Vector{Complex{T}}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real
adams_moulton_outward!(Z::Vector{Complex{T}}, def::Def{T}, adams::Adams{T}) where T<:Real
```

## Radial integration - inward
```@docs
INSCH!(Z::Vector{Complex{T}}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real
INSCH_WKB!(Z::Vector{Complex{T}}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real
adams_moulton_inward!(Z::Vector{Complex{T}}, def::Def{T}, adams::Adams{T}) where T<:Real
```

## Radial integration - boundary condition
```@docs
adams_moulton_normalize!(Z::Vector{Complex{T}}, ΔQ::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real
```

## Adams-Moulton Master procedures
```@docs
adams_moulton_nodes(E::Real, scr::Vector{T}, grid::CamiDiff.Grid{T}, def::Def{T}; imax=25, msg=true) where T<:Real
adams_moulton_iterate!(Z::Vector{Complex{T}}, init::Init{T}, grid::CamiDiff.Grid{T}, def::Def{T}, adams::Adams{T}; imax=25, ϵ=1e-6, msg=true) where T<:Real
adams_moulton_report_nodes(i::Int, init::Init{T}, grid::CamiDiff.Grid{T}, def::Def{T}, strΔT::String; unitIn="Hartree", msg=true) where T<:Real
adams_moulton_report_iterate(i::Int, imax::Int, init::Init{T}, ϵ, grid::CamiDiff.Grid{T}, def::Def{T}, strΔT::String; unitIn="Hartree", msg=true) where T<:Real
```
