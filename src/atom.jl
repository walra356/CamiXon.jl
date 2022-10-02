# ======================== Atom(name, symbol, Z, I, Q, M, I, gI) ===============

"""
    Atom(Z, A, Q, Zc, element, isotope)

Type with fields:
* `      .Z`:  atomic number (`::Int`)
* `      .A`:  atomic mass number in amu (`::Int`)
* `      .Q`:  ionic charge in a.u. (`::Int`)
* `     .Zc`:  Rydberg charge in a.u. (`::Int`)
* `.element`:  (`::Element`)
* `.isotope`:  (`::Isotope`)

The type `Atom` is best created with the function [`castAtom`](@ref).
"""
struct Atom                           # Isotopic properties
    Z::Int             # atomic number
    A::Int             # atomic mass number in amu
    Q::Int             # ionic charge in a.u.
    Zc::Int            # Rydberg charge in a.u.
    element::Element
    isotope::Isotope
end
# ================================ End =========================================

# =========== castAtom(Z, A, Q, msg=true)) ================================

#...............................................................................
function _stdAtom(Z::Int, A::Int, Q::Int)

    dict = dictIsotopes
    atom = (Z,A) ∈ keys(dict) ? castAtom(;Z, A, Q, msg=false) : return nothing

    return atom

end
#...............................................................................
function _strAtom(Z::Int, A::Int, Q::Int)

    dict = dictIsotopes
    atom = (Z,A) ∈ keys(dict) ? castAtom(;Z, A, Q, msg=false) : return nothing

    strQ = abs(Q) > 1 ? sup(abs(Q)) : ""
    strQ = Q > 0 ? (strQ * 'ᐩ') : Q < 0 ? (strQ * 'ᐨ') : ""
    strN = Q ≠ 0 ? " ion" : ", neutral atom"

    str = atom.isotope.name * strN
    str *= ", " * atom.isotope.symbol * strQ
    str *= ", Z=$Z"
    str *= ", A=$A"
    str *= ", Q=$Q"
    str *= ", Zc=$(Q+1)"

    return str

end
#...............................................................................
function _infoAtom(Z::Int, A::Int, Q::Int)

    dict = dictIsotopes
    atom = (Z,A) ∈ keys(dict) ? castAtom(;Z, A, Q, msg=false) : return nothing

    strQ = abs(Q) > 1 ? sup(abs(Q)) : ""
    strQ = Q > 0 ? (strQ * 'ᐩ') : Q < 0 ? (strQ * 'ᐨ') : ""
    strN = Q ≠ 0 ? " ion" : ", neutral atom"

    str = "Atom: " * atom.isotope.name * strN
    str *= "\n  symbol: " * atom.isotope.symbol * strQ
    str *= "\n  atomic charge: Z = $Z"
    str *= "\n  Rydberg charge: Zc = $(Q+1)"

    return println(str)

end
#...............................................................................
"""
    listAtom(Z::Int, A::Int, Q::Int[; fmt=Object])

Properties of atom with atomic number `Z`, atomic mass number `A`,
ionic charge `Q`.

Output options: `fmt` =  `Object` (default), `String`, `Info`.
#### Example:
```
listAtom("H", 3, 0) == listAtom(1, 3, 0)
  true

listAtom(1, 3, 0; fmt=Info)
Element: hydrogen
    symbol: H
    element: tritium
    atomic number: Z = 1
    atomic weight (relative atomic mass): 1.008
```
"""
function listAtom(Z::Int, A::Int, Q::Int; fmt=Object)

    fmt === Object && return _stdAtom(Z, A, Q)
    fmt === String && return _strAtom(Z, A, Q)
    fmt === Info && return _infoAtom(Z, A, Q)

    return error("Error: invalid output type")

end
function listAtom(elt::String, A::Int, Q::Int; fmt=Object)

    dict = dictAtomicNumbers
    Z = (elt) ∈ keys(dict) ? get(dict, elt, nothing) : return nothing

    fmt === Object && return _stdAtom(Z, A, Q)
    fmt === String && return _strAtom(Z, A, Q)
    fmt === Info && return _infoAtom(Z, A, Q)

    return error("Error: invalid output type")

end
#...............................................................................
"""
    listAtoms(Z1::Int, Z2::Int, Q::Int[; fmt=Object])

Properties of atoms with atomic number in the range `Z1:Z3` and
ionic charge `Q`.

Output options: `fmt` =  `Object` (default), `String`, `Info`.
#### Example
```
listAtoms(1,3,0) == listAtoms(1:3,0)
  true

listAtoms(1:1, 0; fmt=Info);
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
            isnothing(next) ? false : push!(o, next)
        end
    end

    return o

end
function listAtoms(itrZ::UnitRange{Int}, Q::Int; fmt=Object)

    return listAtoms(itrZ.start,itrZ.stop, Q; fmt)

end
#...............................................................................

"""
    castAtom(;Z=1, A=1, Q=0, msg=true)
    castAtom(elt::String; A=1, Q=0, msg=true)

Create Atom with fields:
* `      .Z`:  atomic number (`::Int`)
* `      .A`:  atomic mass number in amu (`::Int`)
* `      .Q`:  ionic charge in a.u. (`::Int`)
* `     .Zc`:  Rydberg charge in a.u. (`::Int`)
* `.element`:  (`::Element`)
* `.isotope`:  (`::Isotope`)
#### Examples:
```
castAtom("Rb"; A=87, Q=0, msg=false) == castAtom(Z=37, A=87, Q=0, msg=false)
  true

castAtom(Z=1, A=3, Q=0, msg=false)
  Atom(1, 3, 0, 1, Element("hydrogen", "H", 1.008), Isotope("³T", "tritium",
  1, 3, 2, 1.7591, 3.016049281, 1//2, 1, 12.33, 2.97896246, 0.0, nothing))

atom = castAtom(Z=1, A=3, Q=0, msg=true);
  Element created: H, hydrogen, Z=1, weight=1.008
  Isotope created: ³T, tritium, Z=1, A=3, N=2, R=1.7591, M=3.016049281, I=1/2⁺, μI=2.97896246, Q=0.0, RA=trace, (radioactive)
  Atom created: tritium, neutral atom, ³T, Z=1, A=3, Q=0, Zc=1

atom
  Atom(1, 3, 0, 1, Element("hydrogen", "H", 1.008), Isotope("³T", "tritium",
  1, 3, 2, 1.7591, 3.016049281, 1//2, 1, 12.33, 2.97896246, 0.0, nothing))

atom.isotope.T½
  12.33
```
"""
function castAtom(;Z=1, A=1, Q=0, msg=true)

    element = castElement(;Z, msg)
    isotope = castIsotope(;Z, A, msg)

    msg && println("Atom created: " * listAtom(Z, A, Q; fmt=String) )

    return Atom(Z, A, Q, 1+Q, element, isotope)

end
function castAtom(elt::String; A=1, Q=0, msg=true)

    dict = dictAtomicNumbers
    Z = (elt) ∈ keys(dict) ? get(dict, elt, nothing) :
                return error("Error: element $(elt) - not found in `dictAtomNumbers`")

    element = castElement(;Z, msg)
    isotope = castIsotope(elt; A, msg)

    msg && println("Atom created: " * listAtom(Z, A, Q; fmt=String) )

    return Atom(Z, A, Q, 1+Q, element, isotope)

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
# ================================ End =========================================


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
