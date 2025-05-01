# Coulomb integrals

## Angular integrals

```@docs
a_direct(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)
b_exchange(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)
A_direct(k::Int, l::Int)
B_exchange(k::Int, l::Int, l′::Int)
```

## Radial integrals

#### Coulomb repulsion - direct integral

```@docs
Fk(k::Int, P::Vector{T}, grid::CamiDiff.Grid) where T<:Real
```

#### Coulomb repulsion - exchange integral

```@docs
Gk(k::Int, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid) where T<:Real
```

#### Direct screening potentials

```@docs
UFk(k::Int, P::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
UF(spinorbit1::Spinorbit, spinorbit2::Spinorbit, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
```

#### Exchange screening potentials

```@docs
UGk(k::Int, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
UG(spinorbit1::Spinorbit, spinorbit2::Spinorbit, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
```

## Direct and exchange energies

```@docs
𝒥(spinorbit1::Spinorbit, spinorbit2::Spinorbit, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
𝒦(spinorbit1::Spinorbit, spinorbit2::Spinorbit, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
```