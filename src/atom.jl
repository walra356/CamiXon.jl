# ======================== Atom(name, symbol, Z, I, Q, M, I, gI) ===========

@doc raw"""
    Atom(name, symbol, Z, I, Q, M, I, gI)

Type with fields `name::String` (name of element), `symbol::String` (symbol of element),
`Z::Int` (atomic number), `Q::Int` (ionic charge in a.u.), `M::Float6` (nuclear mass in amu),
`I::Rational{Int}` (nuclear spin), `gI::Float64` (nuclear g-factor)
#### Examples:
```
hydrogen = createAtom("Hydrogen","1H",1,0,1.00782503223,1//2,5.585694713)
symbol = hydrogen.symbol
println("symbol = $(symbol)")
 symbol = 1H

gI = hydrogen.gI
println("nuclear g-factor = $(gI)")
 nuclear g-factor = 5.585694713
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

# ======================== createAtom(name, symbol, Z, I, Q, M, I, gI) ===========

@doc raw"""
    createAtom(name, symbol, Z, Q, M, I, gI)

Create Atom Type with fields `name::String` (name of element), `symbol::String` (symbol of element),
`Z::Int` (atomic number), `Q::Int` (ionic charge in a.u.), `M::Float6` (nuclear mass in amu),
`I::Rational{Int}` (nuclear spin), `gI::Float64` (nuclear g-factor).
#### Examples:
```
createAtom("Hydrogen","1H",1,0,1.00782503223,1//2,5.585694713)
 Atom("Hydrogen", "1H", 1, 0, 1.00782503223, 1//2, 5.585694713)

createAtom("Helium-4","4He",2,0,4.00260325413,1//2,0.0)
 Atom("Helium-4", "4He", 2, 0, 4.00260325413, 1//2, 0.0)

He4II = Atom("Helium-4 ion","4He+",2,1,4.00260325413,1//2,0.0)
 Atom("Helium-4 ion", "4He+", 2, 1, 4.00260325413, 1//2, 0.0)
```
"""
function createAtom(name::String, symbol::String, Z::Int, Q::Int, M::Float64, I::Rational{Int}, gI::Float64)

    return Atom(name, symbol, Z, Q, M, I, gI)

end
