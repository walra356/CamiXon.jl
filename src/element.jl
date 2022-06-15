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
    str *= "\n    symbol: " * element.symbol
    str *= "\n    element: " * element.name
    str *= "\n    atomic number: Z = $Z"
    str *= "\n    atomic weight (relative atomic mass): " * repr(element.weight)

    return println(str)

end
#...............................................................................
"""
    listElement(Z::Int; io=stdout)

Properties of element with atomic number `Z`.

Output options: `stdout` (default), `String`, `Info`.
#### Example:
```
listElement(1; io=Info)
Element: hydrogen
    symbol: H
    element: tritium
    atomic number: Z = 1
    atomic weight (relative atomic mass): 1.008
```
"""
function listElement(Z::Int; io=stdout)

    io === stdout && return _stdElement(Z)
    io === String && return _strElement(Z)
    io === Info && return _infoElement(Z)

    return error("Error: invalid output type")

end
#...............................................................................
"""
    listElements(Z1::Int, Z2::Int; io=stdout)
#### Example
```
listElements(1,3) == listElements(1:3)
  true

listElements(1:3; io=Info)
  Element: hydrogen
    symbol: H
    name: hydrogen
    atomic number: Z = 1
    atomic weight (relative atomic mass): 1.008
  Element: helium
    symbol: He
    name: helium
    atomic number: Z = 2
    atomic weight (relative atomic mass): 4.0026
  Element: lithium
    symbol: Li
    name: lithium
    atomic number: Z = 3
    atomic weight (relative atomic mass): 6.94
```
"""
function listElements(Z1::Int, Z2::Int; io=stdout)

    o = []

    for Z=Z1:Z2
        next = listElement(Z; io)
        isnothing(next) ? false : push!(o, next)
    end

    return o

end
function listElements(itrZ; io=stdout)

    return listElements(itrZ.start,itrZ.stop; io)

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
castElement(;Z=1, msg=true)
  Element created: hydrogen
    symbol: H
    atomic number (Z): 1
    atomic weight (relative atomic mass): 1.008 amu

  Element("hydrogen", "H", 1.008)
```
"""
function castElement(;Z=1, msg=true)

    element = Z ∈ keys(dictElements) ? get(dictElements, Z, nothing) :
              error("Error: element Z = $Z not present in `dictElements`")

    (name, symbol, weight) = element

    msg && println("Element created: " * listElement(Z; io=String) )

    return Element(name, symbol, weight)

end
# ================================ End =========================================
