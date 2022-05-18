# ======================== Element(name, symbol, weight) =======================

"""
    Element(name, symbol, weight)

Type with fields:
* `  .name`:  name of element (`::String`)
* `.symbol`:  symbol of element  (`::String`)
* `.weight`:  relative atomic mass - atomic weight (`::Float64`)

The type `Element` is best created with the function [`castElement`](@ref)).
"""
struct Element           # elemental properties
    name::String         # ionic charge (a.u.)
    symbol::String       # nuclear mass (amu)
    weight::Union{Float64, Nothing}      # relative atomic mass (atomic weight)
end

# ======== Isotope(Z, A, radius, mass, I, π, lifetime, mdm, eqm, ra) ===========

"""
    Isotope(Z, A, radius, mass, I, π, lifetime, mdm, eqm, ra)

Type with fields:
* `     .Z`:  atomic number (`::Int`)
* `     .N`:  neutron number (`::Int`)
* `     .A`:  atomic mass number in amu (`::Int`)
* `     .R`:  rms charge radius in Fermi (`::Float64`)
* `     .M`:  *atomic* mass in amu (`::Float64`)
* `     .I`:  nuclear spin in units of ħ  (`::Rational{Int}`)
* `     .π`:  parity of nuclear state (`::Int`)
* `     .lt`:  lifetime inyears (`::Float64`)
* `     .mdm`: nuclear magnetic dipole moment (`::Float64`)
* `     .eqm`: nuclear electric quadrupole moment (`::Float64`)
* `     .ra`:  relative abundance in % (`::Float64`)

The type `Isotope` is best created with the function [`castIsotope`](@ref)).
"""
struct Isotope              # Isotopic properties
     Z::Int            # atomic number
     N::Int            # neutron number
     A::Int            # atomic mass number (amu)
     R::Float64        # rms charge radius (Fermi)
     M::Float64        # nuclear mass (amu)
     I::Rational{Int}  # nuclear spin in units of ħ
     π::Int            # parity of nuclear state
     lt::Float64       # lifetime (years)
     mdm::Float64      # nuclear magnetic dipole moment
     eqm::Union{Float64, Nothing}       # nuclear electric quadrupole moment
     ra::Union{Float64, Nothing}       # relative abundance (%)
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

function _specsElement(Z::Int, elt)

    (name, symbol, weight) = elt

    strweight = weight ≠ nothing ? "$(weight) amu" : "not available"

    str = "Element created: $(name)
    symbol: $(symbol)
    atomic number: Z = $Z
    atomic weight (relative atomic mass): " * strweight

    return str

end

# =========== castElement(name, symbol, weight) ================================

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

    elt = Z ∈ keys(dictElements) ? get(dictElements, Z, nothing) :
              error("Error: element Z = $Z not present in `dictElements`")

    msg && println(_specsElement(Z, elt) )

    (name, symbol, weight) = elt

    return Element(name, symbol, weight)

end

# ====================== _specsIsotope(Z, A) ===================================

function _specsIsotope(Z::Int, A::Int, iso)              # Isotope properties

    (symbol, radius, mass, I, π, lifetime, mdm, eqm, ra) = iso
    (name, symbol, weight) = get(dictElements, Z, nothing)

    I = typeof(I) ∈ [Float16, Float32, Float64] ? rationalize(I) : I
    lt = lifetime == 1e100 ? "stable" : "$(lifetime) years"

    strRA = ra ≠ nothing ? "$(ra) %" : "trace"
    strmdm = mdm ≠ nothing ? "$(mdm) " : "not available"
    streqm = eqm ≠ nothing ? "$(eqm)" : "not available"

    strIsotope = sup(A) * symbol
    mass = mass * 0.0000010000000000

    name = (Z,A) == (1,1) ? "hydrogen"  :
           (Z,A) == (1,2) ? "deuterium" :
           (Z,A) == (1,3) ? "tritium" : name

    str = "Isotope created: " * strIsotope * "
    element: $(name)
    atomic number: Z = $Z
    neutron number: n = $(A-Z)
    atomic mass number: A =  $A amu
    rms nuclear charge radius: R = $(radius) fm
    atomic mass: mass = $(mass) amu
    nuclear spin: I = $(I) ħ
    parity of nuclear state: π = $π
    lifetime: $(lt)
    nuclear magnetic dipole moment: mdm = " * strmdm * "
    nuclear electric quadrupole moment: eqm = " * streqm * "
    relative abundance: RA = " * strRA

    return str

end

# ================= castIsotope(;Z=1, A=1, msg=true) ===========================

"""
    castIsotope(;Z=1, A=1, msg=true)

Create Isotope with fields
* `     .Z`:  atomic number (`::Int`)
* `     .N`:  neutron number (`::Int`)
* `     .A`:  atomic mass number in amu (`::Int`)
* `     .R`:  rms charge radius in Fermi (`::Float64`)
* `     .M`:  atomic mass in amu (`::Float64`)
* `     .I`:  nuclear spin in units of ħ (`::Rational{Int}`)
* `     .π`:  parity of nuclear state (`::Int`)
* `     .ra`:  relative abundance in % (`::Float64`)
* `     .mdm`: nuclear magnetic dipole moment (`::Float64`)
* `     .eqm`: nuclear electric quadrupole moment (`::Float64`)
* `     .lt`:  lifetime in years (`::Float64`)
#### Examples:
```
isotope = castIsotope(Z=1, A=3, msg=false)
  Isotope(1,-2, 3, 1.7591, 3.016049281, 1//2, 1, 12.33, 2.97896246, 0, nothing)

isotope.ra
  99.9855

castIsotope(Z=1,A=3);
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
```
"""
function castIsotope(;Z=1, A=1, msg=true)

    dict = dictIsotopes
    isotope = (Z, A) ∈ keys(dict) ? get(dict, (Z, A), nothing) :
    error("Error: isotope (Z = $Z, A = $A) not present in `dictIsotopes`")

    msg && println(_specsIsotope(Z, A, isotope) )

    (symbol, radius, mass, I, π, lifetime, mdm, eqm, ra) = isotope

    mass = mass * 0.000001000000

    return Isotope(Z, A-Z, A, radius, mass, I, π, lifetime, mdm, eqm, ra)

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

    msg && println(_specsAtom(Z, A, Q, element) )

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
