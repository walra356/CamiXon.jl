# Coulomb integrals

## Angular integrals

```@docs
a_direct(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)
b_exchange(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)
```

## Radial integrals

#### Direct integral

```@docs
Fk(k::Int, P::Vector{T}, grid::CamiDiff.Grid) where T<:Real
```

#### Exchange integral

```@docs
Gk(k::Int, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid) where T<:Real
```

#### Direct screening potentials

```@docs
UFk(k::Int, P::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
UF(orbit1::Orbit, orbit2::Orbit, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
```

#### Exchange screening potentials

```@docs
UGk(k::Int, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
UG(orbit1::Orbit, orbit2::Orbit, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
```

## Direct and exchange energies

```@docs
𝒥(orbit1::Orbit, orbit2::Orbit, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
𝒦(orbit1::Orbit, orbit2::Orbit, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
```