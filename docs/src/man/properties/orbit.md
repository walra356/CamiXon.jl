# Orbital properties

## Configuration tools
```@docs
extractCore(config::String)
extractValence(config::String)
collectConfig(config::String)
```

## Orbital
```@docs
dictAtomicOrbital
Orbit
castOrbit(;n=1, ℓ=0, mℓ=0, msg=true)
```

## Spinorbital
```@docs
Spinorbit
castSpinorbit(;n=1, ℓ=0, mℓ=0, up=true, restricted=false, msg=true)
collectSpinorbit(config::String; msg=false)
```

## Shells
```@docs
Shell
castShell(;n=1, ℓ=0, msg=false)
Shells
castShells(strShells::String; msg=false)
```

## Term
```@docs
Term
castTerm(n::Int; ℓ=0, S=1//2, L=0, J=1//2, msg=true)
```