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
#                         shell.jl
# ==============================================================================


# ------------------------------------------------------------------------------
#                       Shell(name, N, n, ℓ)
# ------------------------------------------------------------------------------

function _strShell(name, wl, n, ℓ)

    str = "Shell configuration: $(name)
    shell occupation: wl = $(wl) electrons
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
julia> castShell(n=3, ℓ=0, msg=true)
Shell configuration: 3s²
    shell occupation: wl = 2 electrons
    principal quantum number: n = 3
    orbital angular momentum of electrons: ℓ = 0
Shell("3s²", Spinorbit[Spinorbit("3s↓", 3, 2, 0, 0, -1//2), Spinorbit("3s↑", 3, 2, 0, 0, 1//2)])
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

    wl = 2(2ℓ+1)
    strL = ['s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u']
    name = ℓ > 15 ? "[n=$(n), ℓ=$(ℓ)]" : string(n) * strL[ℓ + 1]
    name *= sup(wl)

    msg && println(_strShell(name, wl, n, ℓ) )

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







# ------------------------------------------------------------------------------
#                       Core(name, N, n, ℓ)
# ------------------------------------------------------------------------------

@doc raw"""
    Core(name, shells)

Type for specification of electron core [`Core`](@ref) with fields:
* `.name`: core configuration (`::String`)
* `.count`: number of shells (`::Int`)
* `.n`: array of shell principal quantum numers (`Vector{Int}`)
* `.ℓ`: array of shell angular momenta (`::Vector{Int}`)
* `.shell`: Array of Shells (`::Vector{Shell}`)

The type `Core` is best created with the function `castCore`.
"""
struct Core
    
    name::String
    count::Int
    n::Vector{Int}
    ℓ::Vector{Int}
    shell::Vector{Shell}
    
end

@doc raw"""
    castCore(name::String; msg=false)

Create configuration of closed electron [`Shell`](@ref)s with fields:
* `.name`: core configuration (`::String`)
* `.count`: number of shells (`::Int`)
* `.n`: array of shell principal quantum numers (`Vector{Int}`)
* `.ℓ`: array of shell angular momenta (`::Vector{Int}`)
* `.shell`: Array of Shells (`::Vector{Shell}`)

#### Example:
```

```
"""
function castCore(name::String; msg=false)

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

    return Core(name, k, on, ol, os)
        
end






