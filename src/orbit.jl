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
#                               orbit.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#                       Orbit(name, n, n′, ℓ, up)
# ------------------------------------------------------------------------------

@doc raw"""
    Orbit(name, n, n′, ℓ, mℓ)

Type for specification of *atomic orbitals* with fields:
* `.name`: name
* ` .n`:  principal quantum number
* `.n′`:  radial quantum number (number of nodes in radial wavefunction)
* ` .ℓ`:  orbital angular momentum valence electron
* `.mℓ`:  orbital angular momentum projection valence electron

The type `Orbit` is best created with the function `castOrbit`.
"""
struct Orbit
    name::String          # LS term notation
     n::Int               # principal quantum number
    n′::Int               # radial quantum number (number of nodes)
     ℓ::Int               # orbital angular momentum valence electron
    mℓ::Int               # orbital angular momentum projection valence electron
end

# ------------------------------------------------------------------------------
#                       castOrbit(;n=1, ℓ=0, mℓ=0, msg=true)
# ------------------------------------------------------------------------------

function _specsOrbit(name, n, n′, ℓ, mℓ)

    str = "Orbital: $(name)
    principal quantum number: n = $n
    radial quantum number: n′ = $(n′) (number of nodes in radial wavefunction)
    orbital angular momentum of valence electron: ℓ = $ℓ
    orbital angular momentum projection of valence electron: mℓ = $(mℓ)"

    return str

end

@doc raw"""
    castOrbit(;n=1, ℓ=0, mℓ=0, msg=true)

Create `Orbit` with fields:
* `.name`: name
* ` .n`:  principal quantum number
* `.n′`:  radial quantum number (number of nodes in radial wavefunction)
* ` .ℓ`:  orbital angular momentum valence electron
* `.mℓ`:  orbital angular momentum projection valence electron
#### Example:
```
julia> castOrbit(n=1, ℓ=0; msg=true)
Orbital: 1s
    principal quantum number: n = 1
    radial quantum number: n′ = 0 (number of nodes in radial wavefunction)
    orbital angular momentum of valence electron: ℓ = 0
    orbital angular momentum projection of valence electron: mℓ = 0
Orbit("1s", 1, 0, 0, 0)
```
    castOrbit(strOrbit::String; mℓ=0, msg=false)

#### Example:
```
julia> castOrbit("2p"; mℓ=-1, msg=true)
Orbital: 2p
    principal quantum number: n = 2
    radial quantum number: n′ = 0 (number of nodes in radial wavefunction)
    orbital angular momentum of valence electron: ℓ = 1
    orbital angular momentum projection of valence electron: mℓ = -1
Orbit("2p", 2, 0, 1, -1)
```
"""
function castOrbit(;n=1, ℓ=0, mℓ=0, msg=false)

    ℓ < n || return error("Error: ℓ < n rule not satisfied")
    (-ℓ ≤ mℓ ≤ ℓ) || return error("Error: -ℓ ≤ mℓ ≤ ℓ rule not satisfied")

    strL = ['s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u']

    name = ℓ > 15 ? "[n=$(n), ℓ=$(ℓ)]" : string(n) * strL[ℓ + 1]

    n′ = n - ℓ - 1

    msg && println(_specsOrbit(name, n, n′, ℓ, mℓ) )

    return Orbit(name, n, n′, ℓ, mℓ)

end
function castOrbit(strOrbit::String; mℓ=0, msg=false)

    nl = strip(lowercase(strOrbit))

    n, ℓ = get(dictAtomicOrbital, nl, nothing)

    return castOrbit(;n, ℓ, mℓ, msg)
        
end

# ========== Spinorbit(name::String, orbit::Orbit, ms::Rational(Int)) ===========

@doc raw"""
    Spinorbit

Type for specification of *atomic Spinorbitals* with fields:
* ` .name`: spinorbital name (string)
* ` .n`:  principal quantum number
* `.n′`:  radial quantum number (number of nodes in radial wavefunction)
* ` .ℓ`:  orbital angular momentum valence electron
* `.mℓ`:  orbital angular momentum projection valence electron
* `.ms`: spin magnetic quantum number (Rational{Int})

The type `Spinorbit` is best created with the function `castSpinorbit`.
"""
struct Spinorbit
    name::String         # spinorbital name
    n::Int               # principal quantum number
    n′::Int               # radial quantum number (number of nodes)
     ℓ::Int               # orbital angular momentum valence electron
    mℓ::Int               # orbital angular momentum projection valence electron
    ms::Rational{Int}    # spin magnetic quantum number
end

# ============= castSpinorbit(;n::Int, ℓ::Int, mℓ::Int, up::Bool, msg:Bool) ===========

function _strSpinorbit(name, n, n′, ℓ, mℓ, ms)

    str = "Spinorbital: $(name)
    principal quantum number: n = $n
    radial quantum number: n′ = $(n′) (number of nodes in radial wavefunction)
    orbital angular momentum of valence electron: ℓ = $ℓ
    orbital angular momentum projection of valence electron: mℓ = $(mℓ)
    spin magnetic quantum number: ms = " * strRational(ms)

    return str

end

@doc raw"""
    castSpinorbit(;n=1, ℓ=0, mℓ=0, up=true, msg=false)

Create `Spinorbit` with fields:
* `.name`: spinorbital name (string)
* `   .n`:  principal quantum number
* `  .n′`:  radial quantum number (number of nodes in radial wavefunction)
* `   .ℓ`:  orbital angular momentum valence electron
* `  .mℓ`:  orbital angular momentum projection valence electron
* `  .ms`: spin magnetic quantum number (Rational{Int})
#### Example:
```
julia> castSpinorbit(n=1, ℓ=0, msg=true)
Spinorbital: 1s↑
    principal quantum number: n = 1
    radial quantum number: n′ = 0 (number of nodes in radial wavefunction)
    orbital angular momentum of valence electron: ℓ = 0
    orbital angular momentum projection of valence electron: mℓ = 0
    spin magnetic quantum number: ms = 1/2
Spinorbit("1s↑", 1, 0, 0, 0, 1//2)
```
    castSpinorbit(config::String; msg=false)

#### Example:
```
julia> castSpinorbit1("2p"; msg=true)
Spinorbital: 2p↓
    principal quantum number: n = 2
    radial quantum number: n′ = 0 (number of nodes in radial wavefunction)
    orbital angular momentum of valence electron: ℓ = 1
    orbital angular momentum projection of valence electron: mℓ = 1
    spin magnetic quantum number: ms = -1/2
Spinorbit("2p↓", 2, 0, 1, 1, -1//2)

julia> castSpinorbit1("2p↑-1"; msg=true)
Spinorbital: 2p↑
    principal quantum number: n = 2
    radial quantum number: n′ = 0 (number of nodes in radial wavefunction)
    orbital angular momentum of valence electron: ℓ = 1
    orbital angular momentum projection of valence electron: mℓ = -1
    spin magnetic quantum number: ms = 1/2
Spinorbit("2p↑", 2, 0, 1, -1, 1//2)
```
"""
function castSpinorbit(;n=1, ℓ=0, mℓ=0, ms=1/2, msg=false)

    ℓ < n || return error("Error: ℓ < n rule not satisfied")
    (-ℓ ≤ mℓ ≤ ℓ) || return error("Error: -ℓ ≤ mℓ ≤ ℓ rule not satisfied")

    strL = ['s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u']

    name = ℓ > 15 ? "[n=$(n), ℓ=$(ℓ)]" : string(n) * strL[ℓ + 1]

    n′ = n - ℓ - 1

    ms = rationalize(ms)

    (ms == 1//2) ⊻ (ms == -1//2) || error("Error: unphysical spin (must be 1/2 or -1/2")
    
    name *= string(ms==1/2 ? :↑ : :↓)

    msg && println(_strSpinorbit(name, n, n′, ℓ, mℓ, ms) )

    return Spinorbit(name, n, n′, ℓ, mℓ, ms)

end
function castSpinorbit(config::String; msg=false)

    c = collect(strip(lowercase(config)))
    l = length(c) 
    l < 6 || error("Error: $c not recognized as a spinorbital configuration")
    
    if l == 2
        mℓ = 0
        ms = -1//2
    elseif l == 3
        mℓ = 0
        ms = c[3] == '↑' ? 1//2 : -1//2
    elseif l == 4
        mℓ = parse(Int, c[4]) 
        ms = c[3] == '↑' ? 1//2 : -1//2
    else
        mℓ = -parse(Int, c[5]) 
        ms = c[3] == '↑' ? 1//2 : -1//2
    end

    n, ℓ = get(dictAtomicOrbital, join(c[1:2]), join(c[1:2])*": not recognized ")

    return castSpinorbit(;n, ℓ, mℓ, ms, msg)
        
end

# ======================== Term(name, n, ℓ, S, L, J) ===========

@doc raw"""
    Term(name::String, n::Int, ℓ::Int, S::Real, L::Int, J::Real)

Type for specification of atomic *fine-structure Terms* with fields:
* `name`: name
* ` .n`:  principal quantum number
* `.n′`:  radial quantum number (number of nodes in wavefunction)
* ` .ℓ`:  orbital angular momentum valence electron
* ` .S`:  total electron spin in units of ħ
* ` .L`:  total orbital angular momentum in units of ħ
* ` .J`:  total electronic angular momentum in units of ħ

The type `Term` is best created with the function `castTerm`.
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

@doc raw"""
    castTerm(n::Int; ℓ=0, S=1//2, L=0, J=1//2, msg=false)

Specify Term in the *Term notatation* with fields:
* `.n`: principal quantum number
* `.n′`: radial quantum number (number of nodes - autogenerated)
* `.ℓ`: orbital angular momentum valence electron
* `.S`: total electron spin
* `.L`: total orbital angular momentum
* `.J`: total electronic angular momentum
#### Examples:
```
julia> castTerm(1; ℓ=0, S=1//2, L=0, J=1//2, msg=true)
Term created: 1s ²S₁⸝₂; n = 1,  n′ = 0, ℓ = 0, S = 1//2, L = 0, J = 1//2
Term("1s ²S₁⸝₂", 1, 0, 0, 1//2, 0, 1//2)
```
"""
function castTerm(n::Int; ℓ=0, S=1//2, L=0, J=1//2, msg=false)

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

# ------------------------------------------------------------------------------
#               collectSpinorbit(strCore::String; restricted=true)
# ------------------------------------------------------------------------------

@doc raw"""
    collectSpinorbit(config::String; msg=false)

Collect the spinobitals specified by `config` (in standard configuration notation)
into an array of spinorbitals.
#### Example:
```
julia> config = "1s²2s²2p⁶";

julia> collectSpinorbit(config; msg=true)
10-element Vector{Spinorbit}:
 Spinorbit("1s↓", 1, 0, 0, 0, -1//2)
 Spinorbit("1s↑", 1, 0, 0, 0, 1//2)
 Spinorbit("2s↓", 2, 1, 0, 0, -1//2)
 Spinorbit("2s↑", 2, 1, 0, 0, 1//2)
 Spinorbit("2p↓", 2, 0, 1, -1, -1//2)
 Spinorbit("2p↓", 2, 0, 1, 0, -1//2)
 Spinorbit("2p↓", 2, 0, 1, 1, -1//2)
 Spinorbit("2p↑", 2, 0, 1, -1, 1//2)
 Spinorbit("2p↑", 2, 0, 1, 0, 1//2)
 Spinorbit("2p↑", 2, 0, 1, 1, 1//2)
```
"""
function collectSpinorbit(config::String; msg=false)

    c = collectConfig(config);

    spinorbit = [castSpinorbit(c[i]; msg) for i ∈ eachindex(c)]

    return spinorbit

end




