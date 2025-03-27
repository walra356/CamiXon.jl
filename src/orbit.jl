# SPDX-License-Identifier: MIT

# author: Jook Walraven - 11-11-2024

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

# ------------------------------------------------------------------------------
#                       Shell(name, N, n, ℓ)
# ------------------------------------------------------------------------------

function _strShell(name, N, n, ℓ)

    str = "Shell: $(name)
    number of shell electrons: N = $N
    principal quantum number: n = $n
    orbital angular momentum of electrons: ℓ = $ℓ"

    return str

end

@doc raw"""
    Shell(name, spinorbitals)

Type for specification of *closed electron shells* with fields:
* `.name`: shell configuration (`::String`)
* `.spinorbit`: Array of spinorbitals (`::Vector{Spinorbit}`)

The type `Shell` is best created with the function `castShell`.
"""
struct Shell
    
    name::String
    spinorbit::Vector{Spinorbit}
    
end

@doc raw"""
    castShell(;n=1, ℓ=0, msg=false)

Create closed electron [`Shell`](@ref) with fields:
* `.name`: shell configuration (`::String`)
* `.spinorbit`: Array of Spinorbitals (`::Vector{Spinorbit}`)
#### Example:
```
julia> castShell(n=1, ℓ=0, msg=true)
Shell: 3s²
    number of shell electrons: N = 2
    principal quantum number: n = 3
    orbital angular momentum of electrons: ℓ = 0
Shell("3s²", Spinorbit[Spinorbit("3s↓", Orbit("3s", 3, 2, 0, 0), -1//2), Spinorbit("3s↑", Orbit("3s", 3, 2, 0, 0), 1//2)])
```
    castShell(strShell::String; msg=false)

#### Example:
```
julia> castShell("3s", msg=false)
Shell("3s²", Spinorbit[Spinorbit("3s↓", Orbit("3s", 3, 2, 0, 0), -1//2), Spinorbit("3s↑", Orbit("3s", 3, 2, 0, 0), 1//2)])
```
"""
function castShell(;n=1, ℓ=0, msg=false)

    o = Spinorbit[]
    for ms = -1//2:1//2
        for mℓ=-ℓ:ℓ
            push!(o, castSpinorbit(;n, ℓ, mℓ, ms))
        end
    end

    N = 2(2ℓ+1)
    strL = ['s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u']
    name = ℓ > 15 ? "[n=$(n), ℓ=$(ℓ)]" : string(n) * strL[ℓ + 1]
    name *= sup(N)

    msg && println(_strShell(name, N, n, ℓ) )

    return Shell(name, o)
        
end
function castShell(strShell::String; msg=false)

    nl = strip(lowercase(strShell))[1:2]

    n, ℓ = get(dictAtomicOrbital, nl, nothing)

    return castShell(;n, ℓ, msg)
        
end



# ------------------------------------------------------------------------------
#                       Shells(name, N, n, ℓ)
# ------------------------------------------------------------------------------

@doc raw"""
    Shells(name, shells)

Type for specification of closed electron [`Shells`](@ref) with fields:
* `.name`: shell configuration (`::String`)
* `.count`: number of shells (`::Int`)
* `.n`: array of shell principal quantum numers (`Vector{Int}`)
* `.ℓ`: array of shell angular momenta (`::Vector{Int}`)
* `.shell`: Array of Shells (`::Vector{Shell}`)

The type `Shells` is best created with the function `castShells`.
"""
struct Shells
    
    name::String
    count::Int
    n::Vector{Int}
    ℓ::Vector{Int}
    shell::Vector{Shell}
    
end

@doc raw"""
    castShells(strShells::String; msg=false)

Create configuration of closed electron [`Shells`](@ref) with fields:
* `.name`: shell configuration (`::String`)
* `.count`: number of shells (`::Int`)
* `.n`: array of shell principal quantum numers (`Vector{Int}`)
* `.ℓ`: array of shell angular momenta (`::Vector{Int}`)
* `.shell`: Array of Shells (`::Vector{Shell}`)

#### Example:
```
julia> castShells("1s2s",msg=true);
Shell: 1s²
    number of shell electrons: N = 2
    principal quantum number: n = 1
    orbital angular momentum of electrons: ℓ = 0
Shell: 2s²
    number of shell electrons: N = 2
    principal quantum number: n = 2
    orbital angular momentum of electrons: ℓ = 0
```
"""
function castShells(strShells::String; msg=false)

    nl = ["1s","2s","2p","3s","3p","3d","4s","4p","4d","4f","5s","5p","5d","5f","5g","6s","6p","6d","7s"]   
    #sh = ["[He]", "[Be]", "[Ne]", "[Mg]", "[Ar]", "[Ca]", "[Zn]", "[Kr]", "[Sr]", "[Cd]", "[Xe]", "[Ba]", "[Yb]", "[Hg]", "[Rn]"]
    
 #   for i ∈ eachindex(sh)
 #       if occursin(sh, strShells)
 #           os = Shell[sh[i]]
 #       end
 #   end
    
    name = ""
    os = Shell[]
    on = Int[]
    ol = Int[]
    k = 0
    for i ∈ eachindex(nl)
        if occursin(nl[i], lowercase(strShells))
            n, ℓ = get(dictAtomicOrbital, nl[i], nothing)
            shell = castShell(;n, ℓ, msg)
            push!(os, shell)
            push!(on, ℓ)
            push!(ol, ℓ)
            name *= shell.name
            k +=1
        end
    end

    return Shells(name, k, on, ol, os)
        
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
#                       extractCore(config::String)
# ------------------------------------------------------------------------------

@doc raw"""
    extractCore(config::String)

Extract *core configuration* from `config` given in *compact configuration notation*
#### Example
```
julia> extractCore("[Ar]4s¹")
"1s²2s²2p⁶3s²3p⁶"
```
"""
function extractCore(config::String)

    str = config[1] === '[' ? config[1:4] : return ""

    strCore = get(dictCoreConfiguration, str, "unknown")

    strCore == "unknown" && error("Error: $(str) (unknown 'core bracket')")

    return strCore    
        
end

# ------------------------------------------------------------------------------
#                       extractValence(config::String)
# ------------------------------------------------------------------------------

@doc raw"""
    extractValence(config::String)

Extract *valence configuration* from `config` given in *compact configuration notation*
#### Example
```
julia> extractValence("[Ar]4s¹")
"4s¹"
```
"""
function extractValence(config::String)
    
    strValence = occursin("[", config) ? config[5:end] : config

    return strValence
        
end

# ------------------------------------------------------------------------------
#                collectconfig(config::String)
# ------------------------------------------------------------------------------

@doc raw"""
    collectConfig(config::String)

#### Example:
```
julia> collectConfig("[Be]2p¹") == ["1s↓0", "1s↑0", "2s↓0", "2s↑0", "2p↓-1"]
true

julia> julia> collectConfig("1s↑1s↓02p↑-1") == ["1s↑", "1s↓0", "2p↑-1"]
true
```
"""
function collectConfig(config::String)

    o = String[]

    if occursin('↓', config) | occursin('↑', config)
        c = collect(config)
        i = length(c)
        while i > 0
            j = ((c[i-2] == '↓') ⊻  (c[i-2] == '↑')) ? 3 : 
                ((c[i-1] == '↓') ⊻  (c[i-1] == '↑')) ? 2 : 1
            a = join(c[i-j-1:i])
            i -= j+2
            push!(o, a)
        end
    else
        config = extractCore(config) * extractValence(config)
        c = collect(config)
        i = length(c)
        while i > 0
            N = isascii(c[i-1]) ? undosup(c[i]) : undosup(join(c[i-1:i]))
            i -= isascii(c[i-1]) ? 1 : 2
            nl = join(c[i-1:i])
            n, ℓ = get(dictAtomicOrbital, nl, "unknown")
            i -= 2
            for k = N:-1:1
                ms = k ≤ 2ℓ+1 ? '↓' : '↑'
                ml = k ≤ 2ℓ+1 ? (-ℓ + k - 1) : (-3ℓ + k - 2)
                push!(o, nl*ms*string(ml))
            end
        end
    end

    return reverse(o)

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