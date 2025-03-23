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
#                               atom.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#                  Atom(Z, A, Q, Zc, element, isotope)
# ------------------------------------------------------------------------------

@doc raw"""
    Atom(Z, A, Q, Zc, element, isotope, config)

Type with fields:
* `      .Z`:  atomic number (`::Int`)
* `      .A`:  atomic mass number in amu (`::Int`)
* `      .Q`:  ionic charge in a.u. (`::Int`)
* `     .Zc`:  Rydberg charge in a.u. (`::Int`)
* `.element`:  (`::Element`)
* `.isotope`:  (`::Isotope`)
* ` .config`:  electron configuration (`::String`)

The type `Atom` is best created with the function [`castAtom`](@ref).
"""
struct Atom                           # Isotopic properties
    Z::Int             # atomic number
    A::Int             # atomic mass number in amu
    Q::Int             # ionic charge in a.u.
    Zc::Int            # Rydberg charge in a.u.
    element::Element
    isotope::Isotope
    config::String     # electron configuration
end

# ------------------------------------------------------------------------------
#                       castAtom(;Z=1, A=1, Q=0, msg=false)
# ------------------------------------------------------------------------------

@doc raw"""
    castAtom(;Z=1, A=1, Q=0, msg=false)
    castAtom(elt::String; A=1, Q=0, msg=false)

Create Atom with fields:
* `      .Z`:  atomic number (`::Int`)
* `      .A`:  atomic mass number in amu (`::Int`)
* `      .Q`:  ionic charge in a.u. (`::Int`)
* `     .Zc`:  Rydberg charge in a.u. (`::Int`)
* `.element`:  (`::Element`)
* `.isotope`:  (`::Isotope`)
* ` .config`:  electron configuration (`::String`)

`elt`: symbolic element name
#### Examples:
```
julia> castAtom("Rb"; A=87, Q=0) == castAtom(Z=37, A=87, Q=0)
true

julia> castAtom(Z=1, A=3, Q=0)
Atom(1, 3, 0, 1, Element("hydrogen", "H", 1.008), Isotope("³T", "tritium", 1, 3, 2, 1.7591, 3.016049281, 1//2, 1, 12.33, 2.97896246, 0.0, nothing), "1s¹")

julia> atom = castAtom(Z=1, A=3, Q=0, msg=true);
Element created: H, hydrogen, Z=1, weight=1.008
Isotope created: ³T, tritium, Z=1, A=3, N=2, R=1.7591, M=3.016049281, I=1/2⁺, μI=2.97896246, Q=0.0, RA=trace, (radioactive)
Atom: tritium, neutral atom
  symbol: ³T
  atomic charge: Z = 1
  Rydberg charge: Zc = 1
  electron configuration: 1s¹
Atom created: Atom: tritium, neutral atom
  symbol: ³T
  atomic charge: Z = 1
  Rydberg charge: Zc = 1
  electron configuration: 1s¹

julia> atom
Atom(1, 3, 0, 1, Element("hydrogen", "H", 1.008), Isotope("³T", "tritium", 1, 3, 2, 1.7591, 3.016049281, 1//2, 1, 12.33, 2.97896246, 0.0, nothing), "1s¹")

julia> string(atom.isotope.T½) * " seconds"
"12.33 seconds"
```
"""
function castAtom(;Z=1, A=1, Q=0, msg=false)

    element = castElement(;Z, msg)
    isotope = castIsotope(;Z, A, msg)
     config = get(dictConfiguration, (Z, Q), nothing)
     config = isnothing(config) ? "not provided" : config

    msg && println("Atom created: " * listAtom(Z, A, Q; fmt=Info) )

    return Atom(Z, A, Q, 1+Q, element, isotope, config)

end
function castAtom(elt::String; A=1, Q=0, msg=false)

    dict = dictAtomicNumber
    Z = (elt) ∈ keys(dict) ? get(dict, elt, nothing) :
                return error("Error: element $(elt) - not found in `dictAtomNumbers`")

    return castAtom(;Z, A, Q, msg)

end

# ------------------------------------------------------------------------------
#                       listAtom(Z, A, Q; fmt=Object)
# ------------------------------------------------------------------------------

function _stdAtom(Z::Int, A::Int, Q::Int)

    dict = dictIsotope
    atom = (Z,A) ∈ keys(dict) ? castAtom(;Z, A, Q, msg=false) : return nothing

    return atom

end
#...............................................................................
function _strAtom(Z::Int, A::Int, Q::Int)

    dict = dictIsotope
    atom = (Z,A) ∈ keys(dict) ? castAtom(;Z, A, Q, msg=false) : return nothing

    config = get(dictConfiguration, (Z, Q), nothing)
    config = isnothing(config) ? "not provided" : config

    strQ = abs(Q) > 1 ? sup(abs(Q)) : ""
    strQ = Q > 0 ? (strQ * 'ᐩ') : Q < 0 ? (strQ * 'ᐨ') : ""
    strN = Q ≠ 0 ? " ion" : ", neutral atom"

    str = atom.isotope.name * strN
    str *= ", " * atom.isotope.symbol * strQ
    str *= ", Z=$Z"
    str *= ", A=$A"
    str *= ", Q=$Q"
    str *= ", Zc=$(Q+1)"
    str *= ", config=$(config)"

    return str

end
#...............................................................................
function _infoAtom(Z::Int, A::Int, Q::Int; msg=true)

    dict = dictIsotope
    atom = (Z,A) ∈ keys(dict) ? castAtom(;Z, A, Q, msg=false) : return nothing
  
    config = get(dictConfiguration, (Z, Q), nothing)
    config = isnothing(config) ? "not provided" : config

    strQ = abs(Q) > 1 ? sup(abs(Q)) : ""
    strQ = Q > 0 ? (strQ * 'ᐩ') : Q < 0 ? (strQ * 'ᐨ') : ""
    strN = Q ≠ 0 ? " ion" : ", neutral atom"

    str = "Atom: " * atom.isotope.name * strN
    str *= "\n  symbol: " * atom.isotope.symbol * strQ
    str *= "\n  atomic charge: Z = $Z"
    str *= "\n  Rydberg charge: Zc = $(Q+1)"
    str *= "\n  electron configuration: " * config

    msg && println(str)

    return str

end
#...............................................................................
@doc raw"""
    listAtom(Z::Int, A::Int, Q::Int[; fmt=Object])

Properties of atom with atomic number `Z`, atomic mass number `A`,
ionic charge `Q`.

Output options: `fmt` =  `Object` (default), `String`, `Info`.
#### Example:
```
julia> listAtom("H", 3, 0) == listAtom(1, 3, 0)
true

julia> listAtom(1, 3, 0; fmt=Info)
Atom: tritium, neutral atom
  symbol: ³T
  atomic charge: Z = 1
  Rydberg charge: Zc = 1
```
"""
function listAtom(Z::Int, A::Int, Q::Int; fmt=Object, msg=true)

    return fmt === Object ? _stdAtom(Z, A, Q) :
           fmt === String ? _strAtom(Z, A, Q) :
           fmt === Info ? _infoAtom(Z, A, Q; msg) : throw(DomainError(fmt, "unknown output format"))

end
function listAtom(elt::String, A::Int, Q::Int; fmt=Object, msg=true)

    dict = dictAtomicNumber
    Z = (elt) ∈ keys(dict) ? get(dict, elt, nothing) : return nothing

    return listAtom(Z, A, Q; fmt, msg)

end

# ------------------------------------------------------------------------------
#                       listAtoms(Z1, Z2, Q; fmt=Object)
# ------------------------------------------------------------------------------
@doc raw"""
    listAtoms(Z1::Int, Z2::Int, Q::Int[; fmt=Object])

Properties of atoms with atomic number in the range `Z1:Z3` and
ionic charge `Q`.

Output options: `fmt` =  `Object` (default), `String`, `Info`.
#### Example
```
julia> listAtoms(1,3,0) == listAtoms(1:3,0)
true

julia> listAtoms(1:1, 0; fmt=Info);
Atom: hydrogen, neutral atom
  symbol: ¹H
  atomic charge: Z = 1
  Rydberg charge: Zc = 1
Atom: deuterium, neutral atom
  symbol: ²D
  atomic charge: Z = 1
  Rydberg charge: Zc = 1
Atom: tritium, neutral atom
  symbol: ³T
  atomic charge: Z = 1
  Rydberg charge: Zc = 1
```
"""
function listAtoms(Z1::Int, Z2::Int, Q::Int; fmt=Object)

    o = []

    for Z=Z1:Z2
        for A=1:3Z
            next = listAtom(Z, A, Q; fmt)
            isnothing(next) ? true : Base.push!(o, next)
        end
    end

    return o

end
function listAtoms(itr::UnitRange{Int}, Q::Int; fmt=Object)

    return listAtoms(itr.start, itr.stop, Q; fmt)

end