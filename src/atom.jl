# ======================== Atom(name, symbol, Z, I, Q, M, I, gI) ===========

@doc raw"""
    Atom(name, symbol, Z, I, Q, M, I, gI)

Type to specify the *atomic species* with fields `name::String` (name of element),
`symbol::String` (symbol of element), `Z::Int` (atomic number),
`Q::Int` (ionic charge in a.u.), `M::Float6` (nuclear mass in amu),
`I::Rational{Int}` (nuclear spin), `gI::Float64` (nuclear g-factor)
#### Examples:
```
hydrogen = createAtom(1,0,1.00782503223,1//2,5.585694713)
 Atom("Hydrogen", "¹H", 1, 0, 1.00782503223, 1//2, 5.585694713)

name = hydrogen.name
symbol = hydrogen.symbol
Z = hydrogen.Z
Q = hydrogen.Q
M = hydrogen.M
I = hydrogen.I
g = hydrogen.gI
println("name = $(name), symbol = $(symbol), Z = $Z, Q = $Q, M = $M, I = $I, gI = $g")
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
 Atom("Hydrogen", "¹H", 1, 0, 1.00782503223, 1//2, 5.585694713)

createAtom(2,1,4.00260325413,1//2,0.0)
 Atom("Helium ion", "⁴Heᐩ", 2, 1, 4.00260325413, 1//2, 0.0)
```
"""
function createAtom(Z::Int, Q::Int, M::Float64, I::Rational{Int}, gI::Float64)

    strQ = abs(Q) > 1 ? sup(abs(Q)) : ""
    strQ = Q > 0 ? (strQ * 'ᐩ') : Q < 0 ? (strQ * 'ᐨ') : ""

    (name,symbol) = mendeleev(Z)

    name = Q ≠ 0 ? (name * " ion") : name
    symbol = sup(Int(round(M))) * symbol * strQ

    return Atom(name, symbol, Z, Q, M, I, gI)

end

# ======================== Term(n, ℓ, S, L, J) ===========

@doc raw"""
    Term(n, ℓ, S, L, J)

Type to specify the *fine-structure Term* in *Russell-Saunders notation* with fields
`notation::String`, `n::Int` (principal quantum number),
`ℓ::Int` (orbital angular momentum valence electron),
`S::Rational` (total electron spin), `L::Int` (total orbital angular momentum),
`J::Rational` (total electronic angular momentum).
#### Examples:
```
term_H1I = Term("1s ²S₁⸝₂", 1, 0, 1//2, 0, 1//2)

notation = term_H1I.notation
n = term_H1I.n
ℓ = term_H1I.ℓ
S = term_H1I.S
L = term_H1I.L
J = term_H1I.J
println("notation = $(notation), n = $n, ℓ = $ℓ, S = $S, L = $L, J = $J")
 notation = "1s ²S₁⸝₂", n = 1, ℓ = 0, S = 1//2, L = 0, J = 1//2
```
"""
struct Term
    notation::String     # LS term
    n::Int               # principal quantum number
    ℓ::Int               # orbital angular momentum valence electron
    S::Rational          # total electron spin
    L::Int               # total orbital angular momentum
    J::Rational          # total electronic angular momentum
end

# ======================== createTerm(n, ℓ, S, L, J) ===========

@doc raw"""
    createTerm(n, ℓ, S, L, J)

Specify Term Type of *one* valence electron in the *Term notatation* with fields
`n::Int` (principal quantum number), `ℓ::Int` (orbital angular momentum valence electron),
`S::Rational` (total electron spin), `L::Int` (total orbital angular momentum),
`J::Rational` (total electronic angular momentum).
#### Examples:
```
term_H1I = createTerm(1, 0, 1//2, 0, 1//2)
 Term("1s ²S₁⸝₂", 1, 0, 1//2, 0, 1//2)

notation = term_H1I.notation
n = term_H1I.n
ℓ = term_H1I.ℓ
S = term_H1I.S
L = term_H1I.L
J = term_H1I.J
println("notation = $(notation), n = $n, ℓ = $ℓ, S = $S, L = $L, J = $J")
 notation = "1s ²S₁⸝₂", n = 1, ℓ = 0, S = 1//2, L = 0, J = 1//2
```
"""
function createTerm(n::Int, ℓ::Int, S::Rational, L::Int, J::Rational)

    strL = ['s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u']

    notation = string(n) * strL[ℓ + 1] * ' ' * sup(Int(2S + 1)) * uppercase(strL[L + 1]) * sub(J)

    ℓ < n || return error("jwError: ℓ < n rule not satisfied")
    abs(L-S) ≤ J ≤ (L+S)  || return error("jwError: Δ(LSJ) rule not satisfied")

    return Term(notation, n, ℓ, S, L, J)

end
