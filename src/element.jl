# ======================== Element(name, symbol, weight) =======================

"""
    Element(name, symbol, weight)

Type with fields:
* `  .name`:  name of element (`::String`)
* `.symbol`:  symbol of element  (`::String`)
* `.weight`:  relative atomic mass - atomic weight (`::Float64`)

The type `Element` is best created with the function [`castElement`](@ref).
"""
struct Element           # elemental properties
    name::String         # ionic charge (a.u.)
    symbol::String       # nuclear mass (amu)
    weight::Union{Float64, Nothing}      # relative atomic mass (atomic weight)
end
# ================================ End =========================================

# =========== castElement(name, symbol, weight) ================================

#...............................................................................
function _stdElement(Z::Int)

    dict = dictElements
    element = (Z) ∈ keys(dict) ? castElement(;Z, msg=false) : return nothing

    return element

end
#...............................................................................
function _strElement(Z::Int)

    dict = dictElements
    element = (Z) ∈ keys(dict) ? castElement(;Z, msg=false) : return nothing

    str = element.symbol
    str *= ", " * element.name
    str *= ", Z=$Z"
    str *= ", weight=" * repr(element.weight)

    return str

end
#...............................................................................
function _infoElement(Z::Int)

    dict = dictElements
    element = (Z) ∈ keys(dict) ? castElement(; Z, msg=false) : return nothing

    str = "Element: " * element.name
    str *= "\n  symbol: " * element.symbol
    str *= "\n  atomic number: Z = $Z"
    str *= "\n  atomic weight (relative atomic mass): " * repr(element.weight)

    return println(str)

end
#...............................................................................
"""
    listElement(Z::Int; fmt=Object)

Properties of element with atomic number `Z`.

Output options: `Object` (default), `String`, `Info`.
#### Example:
```
listElement(1; fmt=Info)
  Element: hydrogen
    symbol: H
    atomic number: Z = 1
    atomic weight (relative atomic mass): 1.008
```
"""
function listElement(Z::Int; fmt=Object)

    fmt === Object && return _stdElement(Z, A, Q)
    fmt === String && return _strElement(Z)
    fmt === Info && return _infoElement(Z)

    return error("Error: invalid output type")

end
#...............................................................................
"""
#### Example
```
listElements(1,3) == listElements(1:3)
  true

listElements(1:3; fmt=Info)
  Element: hydrogen
    symbol: H
    atomic number: Z = 1
    atomic weight (relative atomic mass): 1.008
  Element: helium
    symbol: He
    atomic number: Z = 2
    atomic weight (relative atomic mass): 4.0026
  Element: lithium
    symbol: Li
    atomic number: Z = 3
    atomic weight (relative atomic mass): 6.94
```
"""
function listElements(Z1::Int, Z2::Int; fmt=Object)

    o = []

    for Z=Z1:Z2
        next = listElement(Z; fmt)
        isnothing(next) ? false : push!(o, next)
    end

    return o

end
function listElements(itrZ; fmt=Object)

    return listElements(itrZ.start,itrZ.stop; fmt)

end
# ------------------------------------------------------------------------------
"""
    castElement(;Z=1, msg=true)

Create Atom with fields
* `  .name`:  name of element
* `.symbol`:  symbol of element
* `.weight`:  relative atomic mass (atomic weight)
#### Example:
```
element = castElement(;Z=1, msg=true)
element
  Element created: H, hydrogen, Z=1, weight=1.008

  Element("hydrogen", "H", 1.008)
```
"""
function castElement(;Z=1, msg=true)

    element = Z ∈ keys(dictElements) ? get(dictElements, Z, nothing) :
              error("Error: element Z = $Z not present in `dictElements`")

    (name, symbol, weight) = element

    msg && println("Element created: " * listElement(Z; fmt=String) )

    return Element(name, symbol, weight)

end
# ================================ End =========================================
