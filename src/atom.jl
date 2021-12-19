# ======================== Atom(name, symbol, Z, I, Q, M, I, gI) ===============

@doc raw"""
    Atom(name::String, symbol::String, Z::Int, Zc::Int, Q::Int, M::Float64, I::Real, gI::Float64)

Type with fields:
*  `name`:  name of element
*`symbol`:  symbol of element
*     `Z`:  atomic number
*    `Zc`:  Rydberg charge in a.u.
*     `Q`:  ionic charge in a.u.
*     `M`:  nuclear mass in amu
*     `I`:  nuclear spin in units of ħ
*    `gI`:  nuclear g-factor

Note: the type `Atom` is best created by the function `createAtom`.
#### Examples:
```
atom = Atom("Hydrogen", "¹H", 1, 1, 0, 1.00782503223, 1//2, 5.585694713)
atom
 Atom("Hydrogen", "¹H", 1, 1, 0, 1.00782503223, 1//2, 5.585694713)
```
"""
struct Atom              # atomic properties
    name::String         # name of element
    symbol::String       # atomic symbol
    Z::Int               # atomic number
    Zc::Int              # Rydberg charge
    Q::Int               # ionic charge (a.u.)
    M::Float64           # nuclear mass (amu)
    I::Real              # nuclear spin as integer or rational number
    gI::Float64          # nucear g-factor
end

# ======================== createAtom(Z, Q, M, I, gI) ===========

@doc raw"""
    createAtom(Z::Int, Q::Int, M::Float64, I::Real, gI::Float64)

Create Atom Type with fields `name::String` (name of element), `symbol::String` (symbol of element),
`Z::Int` (atomic number), `Zc::Int` (radial quantum number), `Q::Int` (ionic charge in a.u.), `M::Float6` (nuclear mass in amu),
`I::Real` (nuclear spin as integer or rational number), `gI::Float64` (nuclear g-factor).
#### Examples:
```
createAtom(1,0,1.00782503223,1//2,5.585694713)
 Atom created: Hydrogen, symbol = ¹H, Z = 1, Zc = 1, Q = 0, M = 1.00782503223, I = 1//2, gI = 5.585694713
 Atom("Hydrogen", "¹H", 1, 1, 0, 1.00782503223, 1//2, 5.585694713)

createAtom(2,1,4.00260325413,1//2,0.0)
 Atom created: Helium ion, symbol = ⁴Heᐩ, Z = 2, Zc = 3, Q = 1, M = 4.00260325413, I = 1//2, gI = 0.0
 Atom("Helium ion", "⁴Heᐩ", 2, 3, 1, 4.00260325413, 1//2, 0.0)
```
"""
function createAtom(Z::Int, Q::Int, M::Float64, I::Real, gI::Float64)

    S = typeof(I) ∈ [Float16,Float32,Float64] ? rationalize(I) : I

    strQ = abs(Q) > 1 ? sup(abs(Q)) : ""
    strQ = Q > 0 ? (strQ * 'ᐩ') : Q < 0 ? (strQ * 'ᐨ') : ""

    (name,symbol) = mendeleev(Z)

    name = Q ≠ 0 ? (name * " ion") : name
    symbol = sup(Int(round(M))) * symbol * strQ

    Zc = 1 + Q

    println("Atom created: $(name), symbol = $(symbol), Z = $Z, Zc = $(Zc), Q = $Q, M = $M, I = $I, gI = $gI")

    return Atom(name, symbol, Z, Zc, Q, M, I, gI)

end

# ======================== Term(name, n, ℓ, S, L, J) ===========

@doc raw"""
    Term(name::String, n::Int, ℓ::Int, S::Real, L::Int, J::Real)

Type for specification of the atomic *fine-structure Term* with fields:
*`name`: name
*   `n`: principal quantum number
*  `n′`: radial quantum number (number of nodes in wavefunction)
*   `ℓ`: orbital angular momentum valence electron
*   `S`: total electron spin in units of ħ
*   `L`: total orbital angular momentum in units of ħ
*   `J`: total electronic angular momentum in units of ħ

Note: the type `Term` is best created by the function `createTerm`.
#### Examples:
```
Term_H1I = Term("1s ²S₁⸝₂", 1, 0, 0, 1//2, 0, 1//2)
 Term("1s ²S₁⸝₂", 1, 0, 0, 1//2, 0, 1//2)
```
"""
struct Term
    name::String         # LS term notation
    n::Int               # principal quantum number
    n′::Int              # radial quantum number (number of nodes)
    ℓ::Int               # orbital angular momentum valence electron
    S::Real              # total electron spin as integer or rational number
    L::Int               # total orbital angular momentum
    J::Real              # total electronic angular momentum as integer or rational number
end

# ======================== createTerm(n, ℓ, S, L, J) ===========

@doc raw"""
    createTerm(n::Int, ℓ::Int, S::Real, L::Int, J::Real)

Specify Term Type in the *Term notatation* with fields:
`n::Int` (principal quantum number), `ℓ::Int` (orbital angular momentum valence electron),
`S::Rational` (total electron spin), `L::Int` (total orbital angular momentum),
`J::Rational` (total electronic angular momentum).
#### Examples:
```
term_H1I = createTerm(1, 0, 1//2, 0, 1//2)
 Term created: 1s ²S₁⸝₂, n = 1, n′ = 1, ℓ = 0, S = 1//2, L = 0, J = 1//2
 Term("1s ²S₁⸝₂", 1, 0, 0, 1//2, 0, 1//2)
```
"""
function createTerm(n::Int, ℓ::Int, S::Real, L::Int, J::Real)

    S = typeof(S) ∈ [Float16,Float32,Float64] ? rationalize(S) : S
    J = typeof(J) ∈ [Float16,Float32,Float64] ? rationalize(J) : J

    strL = ['s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u']

    name = string(n) * strL[ℓ + 1] * ' ' * sup(Int(2S + 1)) * uppercase(strL[L + 1]) * sub(J)

    ℓ < n || return error("jwError: ℓ < n rule not satisfied")
    abs(L-S) ≤ J ≤ (L+S)  || return error("jwError: Δ(LSJ) condition not satisfied")

    n′ = n - ℓ - 1

    println("Term created: $(name); n = $n,  n′ = $(n′), ℓ = $ℓ, S = $S, L = $L, J = $J")

    return Term(name, n, n′, ℓ, S, L, J)

end

# ======================== bohrformula(atom, term) =============================

@doc raw"""
    bohrformula(Z::Int, n::Int)

Hydrogenic energy (in Hartree a.u.) for *atom* with *atomic number* `Z` and *principal quantum number* `n`.
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
bohrformula(Z::Int, n::Int) = -1/2*(Z/n)^2

# ======================== mendeleev(Z) ========================================

@doc raw"""
    mendeleev(Z::Int)

The properties `name` and `symbol` of the *element* with *atomic number* `Z`.
#### Example:
```
mendeleev(11)
 ("Sodium", "Na")
```
"""
function mendeleev(Z::Int)

    0 < Z < 45 || return error("jwError: element not inplemeted")

    element = Dict(1 => ("Hydrogen", "H"), 2 => ("Helium", "He"),
    3 => ("Lithium", "Li"), 4 => ("Beryllium", "Be"), 5 => ("Boron", "B"), 6 => ("Carbon", "C"),
    7 => ("Oxygen", "O"), 8 => ("Nitrogen", "N"), 9 => ("Fluorine", "F"), 10 => ("Neon", "Ne"),
    11 => ("Sodium", "Na"), 12 => ("Magnesium", "Mg"), 13 => ("Aluminium", "Al"), 14 => ("Silicon", "Si"),
    15 => ("Phosphorus", "P"), 16 => ("Sulphor", "S"), 17 => ("Chlorine", "Cl"), 18 => ("Argon", "Ar"),
    19 => ("Potassium", "K"), 20 => ("Calcium", "Ca"), 21 => ("Scandium", "Sc"), 22 => ("Titanium", "Ti"),
    23 => ("Vanadium", "Va"), 24 => ("Chromium", "Cr"), 25 => ("Manganese", "Cl"),26 => ("Iron", "Fe"),
    27 => ("Cobalt", "Co"), 28 => ("Nickel", "Ni"), 29 => ("Copper", "Cu"), 30 => ("Zinc", "Zn"),
    31 => ("Gallium", "Ga"), 32 => ("Germanium", "Ge"), 33 => ("Arsenic", "As"), 34 => ("Selenium", "Se"),
    35 => ("Bromine", "Br"), 36 => ("Krypton", "Mg"),
    37 => ("Rubidium", "Rb"), 38 => ("Strontium", "Sr"), 39 => ("Yttrium", "Y"), 40 => ("Zirconium", "Zr"),
    41 => ("Niobium", "Nb"), 42 => ("Molybdenium", "Mo"), 43 => ("Technetium", "Tc"), 44 => ("Ruthenium", "Ru"))

    return (name,symbol) = get(element,Z,nothing)

end
