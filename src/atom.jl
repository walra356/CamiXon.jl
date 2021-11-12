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
    Atom(name::String, symbol::String, Z::Int, Q::Int, M::Float64, I::Real, gI::Float64)

Type to specify the *atomic species* with fields `name::String` (name of element),
`symbol::String` (symbol of element), `Z::Int` (atomic number),
`Q::Int` (ionic charge in a.u.), `M::Float6` (nuclear mass in amu),
`I::Rational{Int}` (nuclear spin as integer or rational number), `gI::Float64` (nuclear g-factor)

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
    I::Real              # I: nuclear spin as integer or rational number
    gI::Float64          # gI nucear g-number
end

# ======================== createAtom(Z, Q, M, I, gI) ===========

@doc raw"""
    createAtom(Z::Int, Q::Int, M::Float64, I::Real, gI::Float64)

Create Atom Type with fields `name::String` (name of element), `symbol::String` (symbol of element),
`Z::Int` (atomic number), `Q::Int` (ionic charge in a.u.), `M::Float6` (nuclear mass in amu),
`I::Real` (nuclear spin as integer or rational number), `gI::Float64` (nuclear g-factor).
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
function createAtom(Z::Int, Q::Int, M::Float64, I::Real, gI::Float64)

    S = typeof(I) ∈ [Float16,Float32,Float64] ? rationalize(I) : I

    strQ = abs(Q) > 1 ? sup(abs(Q)) : ""
    strQ = Q > 0 ? (strQ * 'ᐩ') : Q < 0 ? (strQ * 'ᐨ') : ""

    (name,symbol) = mendeleev(Z)

    name = Q ≠ 0 ? (name * " ion") : name
    symbol = sup(Int(round(M))) * symbol * strQ

    println("Atom created: $(name), symbol = $(symbol), Z = $Z, Q = $Q, M = $M, I = $I, gI = $gI")

    return Atom(name, symbol, Z, Q, M, I, gI)

end

# ======================== Term(name, n, ℓ, S, L, J) ===========

@doc raw"""
    Term(name::String, n::Int, ℓ::Int, S::Real, L::Int, J::Real)

Type to specify the *fine-structure Term* in *Russell-Saunders notation* with fields
`name::String`, `n::Int` (principal quantum number),
`ℓ::Int` (orbital angular momentum valence electron),
`S::Real` (total electron spin as integer or rational number),
`L::Int` (total orbital angular momentum),
`J::Real` (total electronic angular momentum as integer or rational number).

Note: the type `Term` is best created by the function `createTerm`.
#### Examples:
```
Term_H1I = Term("1s ²S₁⸝₂", 1, 0, 1//2, 0, 1//2)
 Term("1s ²S₁⸝₂", 1, 0, 1//2, 0, 1//2)

name = Term_H1I.name
n = Term_H1I.n
ℓ = Term_H1I.ℓ
S = Term_H1I.S
L = Term_H1I.L
J = Term_H1I.J
println("name = $(name, n = $n, ℓ = $ℓ, S = $S, L = $L, J = $J")
 name = "1s ²S₁⸝₂", n = 1, ℓ = 0, S = 1//2, L = 0, J = 1//2
```
"""
struct Term
    name::String         # LS term notation
    n::Int               # principal quantum number
    ℓ::Int               # orbital angular momentum valence electron
    S::Real              # total electron spin as integer or rational number
    L::Int               # total orbital angular momentum
    J::Real              # total electronic angular momentum as integer or rational number
end

# ======================== createTerm(n, ℓ, S, L, J) ===========

@doc raw"""
    createTerm(n::Int, ℓ::Int, S::Real, L::Int, J::Real)

Specify Term Type in the *Term notatation* with fields
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
function createTerm(n::Int, ℓ::Int, S::Real, L::Int, J::Real)

    S = typeof(S) ∈ [Float16,Float32,Float64] ? rationalize(S) : S
    J = typeof(J) ∈ [Float16,Float32,Float64] ? rationalize(J) : J

    strL = ['s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u']

    name = string(n) * strL[ℓ + 1] * ' ' * sup(Int(2S + 1)) * uppercase(strL[L + 1]) * sub(J)

    ℓ < n || return error("jwError: ℓ < n rule not satisfied")
    abs(L-S) ≤ J ≤ (L+S)  || return error("jwError: Δ(LSJ) condition not satisfied")

    println("Term created: $(name); n = $n, ℓ = $ℓ, S = $S, L = $L, J = $J")

    return Term(name, n, ℓ, S, L, J)

end

# ======================== Grid(N, h, r0, r, r′) ===============

@doc raw"""
    Grid(N::Int, h::Float64, r0::Float64, r::Vector{Float64}, r′::Vector{Float64})

Type to specify the `Grid` on which the atomic wavefunction is defined, with fields `N::Int` (number of grid points),
`h::Float64` (step size on uniform grid), `r0::Float64` (scale factor for physical grid in a.u.),
`r::Vector{Float64}` (tabulated grid function), r′::Vector{Float64} (tabulated derivative of grid function),
k::Int (Adams-Moulton order), am::Vector{Float}.

Note: the type `Grid` is created by the function `createGrid` which serves to tabulate the grid functions.
"""
struct Grid
    N::Int               # number of grid points
    h::Float64           # stepsize on uniform grid
    r0::Float64          # scale factor for physical grid in a.u.
    r::Vector{Float64}   # tabulated grid function
    r′::Vector{Float64}  # tabulated derivative of grid function
    k::Int               # Adams-Moulton order
    am::Vector{Float}    # Adams-Moulton weight coefficients
end

# ======================== gridfunction(n, h, r0; pmax=6, deriv=0)  ===============

@doc raw"""
    gridfunction(n::Int, h::Float64, r0::Float64; pmax=6, deriv=0)

Transformation from the uniform grid (used for finite-difference calculus) and the nonuniform physical grid,
```math
    f[t,h,r0]=r_0[t+\frac{1}{2!}t^2+\frac{1}{3!}t^3+\frac{1}{4!}t^4+\frac{1}{5!}t^5+\frac{1}{6!}t^6]
```
```math
    f^{′}[t,h,r0]=f[t,h,r_0;deriv=1]=r_0h[1+t+\frac{1}{2!}t^2+\frac{1}{3!}t^3+\frac{1}{4!}t^4+\frac{1}{5!}t^5],
```
where ``t=(n-1)h``.
"""
function gridfunction(n::Int, h::Float64, r0::Float64; pmax=6, deriv=0)
# =================================================================
# my grid function - exponential at small r and algebraic a large r
# gridfunction(n::Int; r0=grid.r0) = r0*(exp((n-1) * h)-1.0) # gridfunction used by Walter Johnson type
# =================================================================

    deriv < pmax || return 0.0

    t = (n-1) * h
    o = deriv > 0 ? 1.0 : 0.0
    v = 1.0

    for p=1:(pmax-deriv)
        v = v * t / p
        o = o + v
    end

    return deriv > 0 ? o * r0 * h : o * r0

end

# ======================== createGrid(N; h=0.01, r0=0.001)   ===============

@doc raw"""
    createGrid(N::Int; h=0.01, r0=0.001, k=8)

Tabulate the gridfunction for an array of `N` points with uniform grid spacing `h` and physical scale factor `r0`.
#### Example:
```
grid = createGrid(3; h=0.01, r0=0.001)
grid.r
3-element Vector{Float64}:
 0.0
 1.0050167084168057e-5
 2.020134002675555e-5
```
"""
function createGrid(N::Int; h=0.01, r0=0.001, k=8)

    r = [gridfunction(n, h, r0) for n=1:N]
    r′= [gridfunction(n, h ,r0; deriv=1) for n=1:N]
    am = create_adams_moulton_weights(k)

    return Grid(N, h, r0, r, r′, k, am)

end



# ============================== matG(n, Etot, atom, grid, scr) ================

@doc raw"""
    matG(n::Int, E::Float64, atom::Atom, grid::Grid, scr::Vector{Float64})

Coupling matrix for a set of two coupled differential equations,
```math
    \frac{dy}{dn}[n]=G[n]\thinspace y[n]\equiv f[n],
```
where ``y=\left(\begin{array}{cc} P\\ Q  \end{array}\right)`` and
``G=\left(\begin{array}{cc} 0 & b[n]\\ c[n] & 0 \end{array}\right)``
"""
function matG(n::Int, E::Float64, atom::Atom, grid::Grid, scr::Vector{Float64})
# ==============================================================================
# matG - coupling matrix - Johnson (2.54)
# ==============================================================================
    G = [0.0 1.0; 0.0 0.0]

    Z = atom.Z
    ℓ = term.ℓ
    r = grid.r[n+1]
    r′= grid.r′[n+1]
    s = scr[n+1]

    a = 0.0
    b = 1.0
    c = fGc(n, r, E, Z, ℓ, s)
    d = 0.0

    G[1,1] = a
    G[1,2] = b * r′
    G[2,1] = c * r′
    G[2,2] = d

    return G

end
