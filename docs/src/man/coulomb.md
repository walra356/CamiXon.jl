# Coulomb integrals

## Angular integrals

```@docs
a_direct(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)
b_exchange(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)
```
## Radial integrals

#### direct integrals

```@docs
UFk(k::Int, P::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
Fk(k::Int, P::Vector{T}, grid::CamiDiff.Grid) where T<:Real
```

#### exchange integrals

```@docs
UGk(k::Int, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
Gk(k::Int, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid) where T<:Real
```

## Screening potentials

```@docs
UF(orbit1::Orbit, orbit2::Orbit, P::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
UG(orbit1::Orbit, orbit2::Orbit, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
```

## Direct and exchange integrals

```@docs
𝒥(orbit1::Orbit, orbit2::Orbit, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
𝒦(orbit1::Orbit, orbit2::Orbit, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
```