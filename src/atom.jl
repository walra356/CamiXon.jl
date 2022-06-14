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

# ======================== Atom(name, symbol, Z, I, Q, M, I, gI) ===============

"""
    Atom(Z, A, Q, Zc, element, isotope)

Type with fields:
* `     .Z`:  atomic number (`::Int`)
* `     .A`:  atomic mass number in amu (`::Int`)
* `     .Q`:  ionic charge in a.u. (`::Int`)
* `    .Zc`:  Rydberg charge in a.u. (`::Int`)
* `.element`:  (`::Element`)
* `.isotope`:  (`::Isotope`)

The type `Atom` is best created with the function [`castAtom`](@ref)).
"""
struct Atom                           # Isotopic properties
    Z::Int             # atomic number
    A::Int             # atomic mass number in amu
    Q::Int             # ionic charge in a.u.
    Zc::Int            # Rydberg charge in a.u.
    element::Element
    isotope::Isotope
end

# ======================== Orbit(name, n, n′, ℓ, up) ===========================

"""
    Orbit(name, n, n′, ℓ)

Type for specification of *atomic orbitals* with fields:
* `.name`: name
* ` .n`:  principal quantum number
* `.n′`:  radial quantum number (number of nodes in radial wavefunction)
* ` .ℓ`:  orbital angular momentum valence electron

The type `Orbit` is best created with the function `castOrbit`.
"""
struct Orbit
    name::String         # LS term notation
    n::Int               # principal quantum number
    n′::Int              # radial quantum number (number of nodes)
    ℓ::Int               # orbital angular momentum valence electron
end

# ============================ _specsElement(Z) ================================

function _infoElement(Z::Int, elt)

    (name, symbol, weight) = elt

    strweight = weight ≠ nothing ? "$(weight) amu" : "not available"

    str = "Element created: $(name)
    symbol: $(symbol)
    atomic number: Z = $Z
    atomic weight (relative atomic mass): " * strweight

    return str

end

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
#...............................................................................
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

# ======================= castAtom(;Z=1, A=1, Q=0, msg=true) ===================

function _specsAtom(Z::Int, A::Int, Q::Int)

    (name, symbol, weight) = get(dictElements, Z, nothing)

    name = (Z,A) == (1,2) ? "deuterium" :
           (Z,A) == (1,3) ? "tritium"   : name

    strQ = abs(Q) > 1 ? sup(abs(Q)) : ""
    strQ = Q > 0 ? (strQ * 'ᐩ') : Q < 0 ? (strQ * 'ᐨ') : ""
    strN = Q ≠ 0 ? " ion" : ", neutral atom"

    strD = "D ≡ " * sup(A) * symbol * strQ
    strT = "T ≡ " * sup(A) * symbol * strQ

    symbol = (Z,A) == (1,2) ? strD :
             (Z,A) == (1,3) ? strT : sup(A) * symbol * strQ

    name = name * strN

    charge = Q ≠ 0 ? "ionic charge: $Q" : "neutral atom"

    str = "Atom created: $(name)
    symbol: $(symbol)
    ionic charge: Q = $Q
    Rydberg charge: Zc = $(1+Q)"

    return str

end

"""
    castAtom(;Z=1, A=1, Q=0, msg=true)

Create Atom with fields:
* `     .Z`:  atomic number (`::Int`)
* `     .A`:  atomic mass number in amu (`::Int`)
* `     .Q`:  ionic charge in a.u. (`::Int`)
* `    .Zc`:  Rydberg charge in a.u. (`::Int`)
*`.element`:  (`::Element`)
*`.isotope`:  (`::Isotope`)
#### Examples:
```
atom = castAtom(Z=1, A=1, Q=0, msg=false)
  Atom(1, 1, 0, 1, Element("hydrogen", "H", 1.008), Isotope(1, 0, 1, 0.8783,
  1.007825032, 1//2, 1, 1.0e100, 2.792847351, 0.0, 99.9855))

atom.isotope.ra
  99.9855

castAtom(Z=1, A=3, Q=0);
  Element created: hydrogen
      symbol: H
      atomic number (Z): 1
      atomic weight (relative atomic mass): 1.008 amu
  Isotope created: ³H
      element: tritium
      atomic number: Z = 1
      neutron number: n = 2
      atomic mass number: A =  3 amu
      rms nuclear charge radius: R = 1.7591 fm
      atomic mass: mass = 3.016049281 amu
      nuclear spin: I = 1//2 ħ
      parity of nuclear state: π = 1
      lifetime: 12.33 years
      nuclear magnetic dipole moment: mdm = 2.97896246
      nuclear electric quadrupole moment: eqm = 0
      relative abundance: RA = trace
Atom created: hydrogen - ³H (Z = 1, Zc = 1, Q = 0)
```
"""
function castAtom(;Z=1, A=1, Q=0, msg=true)

    element = castElement(;Z, msg)
    isotope = castIsotope(;Z, A, msg)

    msg && println(_specsAtom(Z, A, Q) )

    return Atom(Z, A, Q, 1+Q, element, isotope)

end

# ======================== castOrbital(n::Int, ℓ::Int) ===========

function _specsOrbit(name, n, n′, ℓ)

    str = "Orbital: $(name)
    principal quantum number: n = $n
    radial quantum number: n′ = $(n′) (number of nodes in radial wavefunction)
    orbital angular momentum of valence electron: ℓ = $ℓ"

    return str

end

"""
    castOrbit(;n=1, ℓ=0, msg=true)

Create `Orbit` with fields:
* `.name`: name
* ` .n`:  principal quantum number
* `.n′`:  radial quantum number (number of nodes in radial wavefunction)
* ` .ℓ`:  orbital angular momentum valence electron
#### Examples:
```
castOrbit(n=1, ℓ=0)
 Orbit created: 1s (n = 1, n′ = 0, ℓ = 0)
 Orbit("1s", 1, 0, 0)
```
"""
function castOrbit(;n=1, ℓ=0, msg=true)

    ℓ < n || return error("Error: ℓ < n rule not satisfied")

    strL = ['s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u']

    name = ℓ > 15 ? "[n=$(n), ℓ=$(ℓ)]" : string(n) * strL[ℓ + 1]

    n′ = n - ℓ - 1

    msg && println(_specsOrbit(name, n, n′, ℓ) )

    return Orbit(name, n, n′, ℓ)

end

# ======================== SpinOrbit(name, n, n′, ℓ, ms) ===========

"""
    SpinOrbit

Type for specification of *atomic spinorbitals* with fields:
* `.name`: name
* ` .n`:  principal quantum number
* `.n′`:  radial quantum number (number of nodes in radial wavefunction)
* ` .ℓ`:  orbital angular momentum valence electron
* `.ms`:  spin magnetic quantum number

The type `SpinOrbit` is best created with the function `createSpinOrbit`.
"""
struct SpinOrbit
    name::String         # LS term notation
    n::Int               # principal quantum number
    n′::Int              # radial quantum number (number of nodes)
    ℓ::Int               # orbital angular momentum valence electron
    ms::Rational{Int}    # spin magnetic quantum number
end


# ====================== createSpinOrbit(o::Orbital; up=true) ==================

"""
    createSpinOrbital(o::Orbit; up=true, msg=true)

Specify `SpinOrbit` with fields:
* `.name`: name
* `   .n`: principal quantum number
* `  .n′`: radial quantum number (number of nodes in radial wavefunction)
* `   .ℓ`: orbital angular momentum valence electron
* `  .ms`: spin magnetic quantum number
#### Examples:
```
s1s = castOrbit(1,0)
createSpinOrbit(s1s; up=true)
  SpinOrbit created: 1s↑ (n = 1, n′ = 0, ℓ = 0, ms = 1//2)
  SpinOrbit("1s↑", 1, 0, 0, 1//2)
```
"""
function createSpinOrbit(o::Orbit; up=true, msg=true)

    name = o.name * string(up ? :↑ : :↓)

    msg && println("SpinOrbit created: $(name) (n = $(o.n), n′ = $(o.n′), ℓ = $(o.ℓ), ms = $(up ? 1//2 : -1//2))")

    return SpinOrbit(name, o.n, o.n′, o.ℓ, (up ? 1//2 : -1//2))

end


# ======================== Term(name, n, ℓ, S, L, J) ===========

"""
    Term(name::String, n::Int, ℓ::Int, S::Real, L::Int, J::Real)

Type for specification of atomic *fine-structure Terms* with fields:
* `name`: name
* ` .n`:  principal quantum number
* `.n′`:  radial quantum number (number of nodes in wavefunction)
* ` .ℓ`:  orbital angular momentum valence electron
* ` .S`:  total electron spin in units of ħ
* ` .L`:  total orbital angular momentum in units of ħ
* ` .J`:  total electronic angular momentum in units of ħ

The type `Term` is best created with the function `createTerm`.
"""
struct Term
    name::String         # LS term notation
    n::Int               # principal quantum number
    n′::Int              # radial quantum number (number of nodes)
    ℓ::Int               # orbital angular momentum valence electron
    S::Real              # total electron spin as integer or rational number
    L::Int               # total orbital angular momentum
    J::Real              # total electronic angular momentum
end

# ===================== createTerm(n::Int; ℓ=0, S=1//2, L=0, J=1//2) ===========

"""
    createTerm(n::Int; ℓ=0, S=1//2, L=0, J=1//2, msg=true)

Specify Term in the *Term notatation* with fields:
* `.n`: principal quantum number
* `.n′`: radial quantum number (number of nodes - autogenerated)
* `.ℓ`: orbital angular momentum valence electron
* `.S`: total electron spin
* `.L`: total orbital angular momentum
* `.J`: total electronic angular momentum
#### Examples:
```
term_H1I = createTerm(1; ℓ=0, S=1//2, L=0, J=1//2)
 Term created: 1s ²S₁⸝₂, n = 1, n′ = 0, ℓ = 0, S = 1//2, L = 0, J = 1//2
 Term("1s ²S₁⸝₂", 1, 0, 0, 1//2, 0, 1//2)
```
"""
function createTerm(n::Int; ℓ=0, S=1//2, L=0, J=1//2, msg=true)

    S = typeof(S) ∈ [Float16,Float32,Float64] ? rationalize(S) : S
    J = typeof(J) ∈ [Float16,Float32,Float64] ? rationalize(J) : J

    strL = ['s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u']

    name = string(n) * strL[ℓ + 1] * ' ' * sup(Int(2S + 1)) * uppercase(strL[L + 1]) * sub(J)

    ℓ < n || return error("Error: ℓ < n rule not satisfied")
    abs(L-S) ≤ J ≤ (L+S)  || return error("Error: Δ(LSJ) condition not satisfied")

    n′ = n - ℓ - 1

    msg && println("Term created: $(name); n = $n,  n′ = $(n′), ℓ = $ℓ, S = $S, L = $L, J = $J")

    return Term(name, n, n′, ℓ, S, L, J)

end

# ======================== bohrformula(atom, term) =============================

@doc raw"""
    bohrformula(Z::Int, n::Int)

Hydrogenic energy (in Hartree a.u.) for *atom* with *atomic number* `Z` and
*principal quantum number* `n`.
```math
    E_n = - \frac{Z^2}{2n^2}
```
#### Example:
```
Z = 2
n = 4
bohrformula(Z,n)
 -0.125
```
"""
bohrformula(Z::Int, n::Int) = -(1//2)*(Z//n)^2
