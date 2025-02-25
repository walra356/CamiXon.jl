# Principal properties

## Element
```@docs
Element
castElement(;Z=1, msg=true)
listElement(Z::Int; fmt=Object)
listElements(Z1::Int, Z2::Int; fmt=Object)
```

## Isotope
```@docs
Isotope
castIsotope(;Z=1, A=1, msg=true)
listIsotope(Z::Int, A::Int; fmt=Object, msg=true)
listIsotopes(Z1::Int, Z2::Int; fmt=Object)
latexIsotopeTable(Z1::Int, Z2::Int; continuation=false)
```

## Atom
```@docs
Atom
castAtom(;Z=1, A=1, Q=0, msg=true)
listAtom(Z::Int, A::Int, Q::Int; fmt=Object)
listAtoms(Z1::Int, Z2::Int, Q::Int; fmt=Object)
```

## Thermodynamic properties
```@docs
melting_point(atomicnumber::Int)
svp(atomicnumber::Int, temp::Real)
latent_heat_vaporization(atomicnumber::Int, temp::Real)
```