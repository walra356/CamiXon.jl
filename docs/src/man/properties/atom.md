# Principal properties

## Element

```@docs
Element
castElement(;Z=1, msg=true)
listElement(Z::Int; fmt=Object)
listElements(Z1::Int, Z2::Int; fmt=Object)
```
#### Dictionary
```@docs
dictElement
```

## Isotope
```@docs
Isotope
castIsotope(;Z=1, A=1, msg=true)
listIsotope(Z::Int, A::Int; fmt=Object, msg=true)
listIsotopes(Z1::Int, Z2::Int; fmt=Object)
latexIsotopeTable(Z1::Int, Z2::Int; continuation=false)
```
#### Dictionary
```@docs
dictIsotope
```

## Atom

```@docs
Atom
castAtom(;Z=1, A=1, Q=0, msg=true)
listAtom(Z::Int, A::Int, Q::Int; fmt=Object)
listAtoms(Z1::Int, Z2::Int, Q::Int; fmt=Object)
```
#### Dictionaries
```@docs
dictAtomicNumber
dictClosedShell
dictConfiguration
```