# Hydrogen

## Atomic properties

```@docs
bohrformula(Z::Int, n::Int)
hydrogenic_reduced_wavefunction(atom::Atom, orbit::Orbit, grid::CamiDiff.Grid{T}) where T<:Real
reduce_wavefunction(Z::Vector{Complex{T}}, grid::CamiDiff.Grid{T}) where T<:Real
restore_wavefunction(Z::Vector{Complex{T}}, atom::Atom, orbit::Orbit, grid::CamiDiff.Grid{T}) where T<:Real
```
## Some special cases
```@docs
RH1s(Z::Int, r::T) where T<:Real
RH2s(Z::Int, r::T) where T<:Real
RH2p(Z::Int, r::T) where T<:Real
```

## Molecular properties
```@docs
silvera_goldman_triplet(r::T) where T<:Real
silvera_goldman_singlet(r::T) where T<:Real
silvera_goldman_exchange(r::T) where T<:Real
silvera_goldman_potential(grid::CamiDiff.Grid{T}; S=0) where T<:Real
rotbarrier(grid::CamiDiff.Grid{T}; ℓ=0) where T<:Real
```