# SPDX-License-Identifier: MIT

# Copyright (c) 2025 Jook Walraven <69215586+walra356@users.noreply.github.com> and contributors

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# ==============================================================================
#                               element.jl
# ==============================================================================

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

    dict = dictElement
    element = (Z) ∈ keys(dict) ? castElement(;Z, msg=false) : return nothing

    return element

end
#...............................................................................
function _strElement(Z::Int)

    dict = dictElement
    element = (Z) ∈ keys(dict) ? castElement(;Z, msg=false) : return nothing

    str = element.symbol
    str *= ", " * element.name
    str *= ", Z=$Z"
    str *= ", weight=" * repr(element.weight)

    return str

end
#...............................................................................
function _infoElement(Z::Int; msg=true)

    dict = dictElement
    element = (Z) ∈ keys(dict) ? castElement(; Z, msg=false) : return nothing

    str = "Element: " * element.name
    str *= "\n  symbol: " * element.symbol
    str *= "\n  atomic number: Z = $Z"
    str *= "\n  atomic weight (relative atomic mass): " * repr(element.weight)

    return msg ? println(str) : str

end
#...............................................................................
"""
    listElement(Z::Int[; fmt=Object])
    listElement(elt::String[; fmt=Object])

Properties of element with atomic number `Z` and symbolic name `elt`

Output options: `fmt` =  `Object` (default), `String`, `Info`.
#### Example:
```
julia> listElement("H") == listElement(1)
true

julia> listElement(1; fmt=Info)
Element: hydrogen
  symbol: H
  atomic number: Z = 1
  atomic weight (relative atomic mass): 1.008
```
"""
function listElement(Z::Int; fmt=Object, msg=true)

    strErr = "Error: invalid output type"

    return fmt === Object ? _stdElement(Z) : fmt === String ? _strElement(Z) :
           fmt === Info ? _infoElement(Z; msg) : error(strErr)

end
function listElement(elt::String; fmt=Object, msg=true)

    dict = dictAtomicNumber
    Z = (elt) ∈ keys(dict) ? get(dict, elt, nothing) : return nothing

    return listElement(Z; fmt,msg)

end
#...............................................................................
"""
    listElements(Z1::Int, Z2::Int[; fmt=Object])
    listElements(itr::UnitRange{Int}; fmt=Object)

Properties of elements with atomic number in the range `itr = Z1:Z2`.

Output options: `fmt` =  `Object` (default), `String`, `Info`.

#### Example
```
julia> listElements(1,3) == listElements(1:3)
true

julia> listElements(1:3; fmt=Info);
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

julia> listElements(1:3; fmt=String)
3-element Vector{Any}:
 "H, hydrogen, Z=1, weight=1.008"
 "He, helium, Z=2, weight=4.0026"
 "Li, lithium, Z=3, weight=6.94"    
```
"""
function listElements(Z1::Int, Z2::Int; fmt=Object)

    o = []

    for Z=Z1:Z2
        next = listElement(Z; fmt)
        isnothing(next) ? true : Base.push!(o, next)
    end

    return o

end
function listElements(itr::UnitRange{Int}; fmt=Object)

    return listElements(itr.start,itr.stop; fmt)

end
# ------------------------------------------------------------------------------
"""
    castElement(;Z=1, msg=true)
    castElement(elt::String; msg=true)

Create Atom with fields
* `  .name`:  name of element
* `.symbol`:  symbol of element
* `.weight`:  relative atomic mass (atomic weight)

  `Z`: atomic number (nuclear charge number)
`elt`: symbolic element name
#### Example:
```
julia> castElement("Rb"; msg=false) == castElement(Z=37, msg=false)
true

julia> element = castElement(;Z=1, msg=true);
Element created: H, hydrogen, Z=1, weight=1.008

julia> element = castElement(;Z=1, msg=false)
Element("hydrogen", "H", 1.008)
```
"""
function castElement(;Z=1, msg=true)

    element = Z ∈ keys(dictElement) ? get(dictElement, Z, nothing) :
              error("Error: element Z = $Z not present in `dictElement`")

    (name, symbol, weight) = element

    msg && println("Element created: " * listElement(Z; fmt=String) )

    return Element(name, symbol, weight)

end
function castElement(elt::String; msg=true)

    dict = dictAtomicNumber
    Z = (elt) ∈ keys(dict) ? get(dict, elt, nothing) :
                return error("Error: element $(elt) - not found in `dictAtomicNumber`")

    dict = dictElement
    element = Z ∈ keys(dict) ? get(dict, Z, nothing) :
              return error("Error: element Z = $Z - not found in `dictElement`")

    (name, symbol, weight) = element

    msg && println("Element created: " * listElement(Z; fmt=String) )

    return Element(name, symbol, weight)

end

# ================================ End =========================================
