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

# ======================== Atom(name, symbol, Z, I, Q, M, I, gI) ===============

@doc raw"""
    Atom(name, symbol, Z, I, Q, M, I, gI)

Type to specify the *atomic species* with fields `name::String` (name of element),
`symbol::String` (symbol of element), `Z::Int` (atomic number),
`Q::Int` (ionic charge in a.u.), `M::Float6` (nuclear mass in amu),
`I::Rational{Int}` (nuclear spin), `gI::Float64` (nuclear g-factor)

Note: the type `Atom` is best created by the function `createAtom`.
#### Examples:
```
Hydrogen = Atom("Hydrogen", "¹H", 1, 0, 1.00782503223, 1//2, 5.585694713)

name = Hydrogen.name
symbol = Hydrogen.symbol
Z = Hydrogen.Z
Q = Hydrogen.Q
M = Hydrogen.M
I = Hydrogen.I
gI = Hydrogen.gI
println("name = $(name), symbol = $(symbol), Z = $Z, Q = $Q, M = $M, I = $I, gI = $gI")
 name = Hydrogen, symbol = ¹H, Z = 1, Q = 0, M = 1.00782503223, I = 1//2, gI = 5.585694713
```
"""
struct Atom              # atomic properties
    name::String         # name of element
    symbol::String       # atomic symbol
    Z::Int               # Z: atomic number
    Q::Int               # Q: ionic charge (a.u.)
    M::Float64           # M: nuclear mass (amu)
    I::Rational{Int}     # I: nuclear spin
    gI::Float64          # gI nucear g-number
end

# ======================== createAtom(Z, Q, M, I, gI) ===========

@doc raw"""
    createAtom(Z, Q, M, I, gI)

Create Atom Type with fields `name::String` (name of element), `symbol::String` (symbol of element),
`Z::Int` (atomic number), `Q::Int` (ionic charge in a.u.), `M::Float6` (nuclear mass in amu),
`I::Rational{Int}` (nuclear spin), `gI::Float64` (nuclear g-factor).
#### Examples:
```
createAtom(1,0,1.00782503223,1//2,5.585694713)
 Atom created: Hydrogen, symbol = ¹H, Z = 1, Q = 0, M = 1.00782503223, I = 1//2, gI = 5.585694713
 Atom("Hydrogen", "¹H", 1, 0, 1.00782503223, 1//2, 5.585694713)

createAtom(2,1,4.00260325413,1//2,0.0)
 Atom created: Helium ion, symbol = ⁴Heᐩ, Z = 2, Q = 1, M = 4.00260325413, I = 1//2, gI = 0.0
 Atom("Helium ion", "⁴Heᐩ", 2, 1, 4.00260325413, 1//2, 0.0)
```
"""

# ======================== Term(note, n, ℓ, S, L, J) ===========

@doc raw"""
    Term(note::String, n::Int, ℓ::Int, S::Rational, L::Int, J::Rational)

Type to specify the *fine-structure Term* in *Russell-Saunders notation* with fields
`note::String`, `n::Int` (principal quantum number),
`ℓ::Int` (orbital angular momentum valence electron),
`S::Rational` (total electron spin), `L::Int` (total orbital angular momentum),
`J::Rational` (total electronic angular momentum).

Note: the type `Term` is best created by the function `createTerm`.
#### Examples:
```
Term_H1I = Term("1s ²S₁⸝₂", 1, 0, 1//2, 0, 1//2)
 Term("1s ²S₁⸝₂", 1, 0, 1//2, 0, 1//2)

note = Term_H1I.note
n = Term_H1I.n
ℓ = Term_H1I.ℓ
S = Term_H1I.S
L = Term_H1I.L
J = Term_H1I.J
println("note = $(note, n = $n, ℓ = $ℓ, S = $S, L = $L, J = $J")
 note = "1s ²S₁⸝₂", n = 1, ℓ = 0, S = 1//2, L = 0, J = 1//2
```
"""
struct Term
    note::String         # LS term notation
    n::Int               # principal quantum number
    ℓ::Int               # orbital angular momentum valence electron
    S::Rational          # total electron spin
    L::Int               # total orbital angular momentum
    J::Rational          # total electronic angular momentum
end

# ======================== createTerm(n, ℓ, S, L, J) ===========

@doc raw"""
    createTerm(n::Int, ℓ::Int, S::Rational, L::Int, J::Rational)

Specify Term Type of *one* valence electron in the *Term notatation* with fields
`n::Int` (principal quantum number), `ℓ::Int` (orbital angular momentum valence electron),
`S::Rational` (total electron spin), `L::Int` (total orbital angular momentum),
`J::Rational` (total electronic angular momentum).
#### Examples:
```
term_H1I = createTerm(1, 0, 1//2, 0, 1//2)
 Term created: 1s ²S₁⸝₂, n = 1, ℓ = 0, S = 1//2, L = 0, J = 1//2
 Term("1s ²S₁⸝₂", 1, 0, 1//2, 0, 1//2)
```
"""
function createTerm(n::Int, ℓ::Int, S::Rational, L::Int, J::Rational)

    strL = ['s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u']

    note = string(n) * strL[ℓ + 1] * ' ' * sup(Int(2S + 1)) * uppercase(strL[L + 1]) * sub(J)

    ℓ < n || return error("jwError: ℓ < n rule not satisfied")
    abs(L-S) ≤ J ≤ (L+S)  || return error("jwError: Δ(LSJ) condition not satisfied")

    println("Term created: $(note); n = $n, ℓ = $ℓ, S = $S, L = $L, J = $J")

    return Term(note, n, ℓ, S, L, J)

end
