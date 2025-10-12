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
#                               config.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#                       formatConfig(config)
# ------------------------------------------------------------------------------

@doc raw"""
    formatConfig(config::String)

format configuration in hydrogenic ordering
#### Example
```
julia> formatConfig("[Ar]4s¹")
"1s²2s²2p⁶3s²3p⁶"
```
"""
function formatConfig(config::String; Err=true)

    o1, o2 = splitConfig(config; Err=Err)
    o = join(o1) * " " * join(o2)

    return o

end

# ------------------------------------------------------------------------------
#                       expandConfig(config; Err=true)
# ------------------------------------------------------------------------------

@doc raw"""
    expandConfig(config::String; Err=true)

Expand configuration into Array of shell configurations.
#### Example
```
julia> expandConfig("[Ar]4s¹") == ["1s²", "2s²", "2p⁶", "3s²", "3p⁶", "4s¹"]
true
```
"""
function expandConfig(config::String; Err=true)
    
        o1, o2 = splitConfig(config; Err=Err)
    
        o = append!(o1, o2)
    
        return o
end




# ------------------------------------------------------------------------------
#                       splitConfig(config; Err=true)
# ------------------------------------------------------------------------------

@doc raw"""
    splitConfig(config::String; Err=true)

Expand configuration into Array of shell configurations.
#### Example
```
julia> splitConfig("[Ar]4s¹") == (["1s²", "2s²", "2p⁶", "3s²", "3p⁶"], ["4s¹"])
true
```
"""
function splitConfig(config::String; Err=true)

    nl = ["1s","2s","2p","3s","3p","3d","4s","4p","4d","4f","5s","5p","5d","5f","5g","6s","6p","6d","7s"]   

    if config[1] === '[' 
        core = get(dictCoreConfiguration, config[1:4], "unknown")
        Err && core == "unknown" && error("Error: $(config[1:4]) (unknown 'core bracket')")
        config = core * config[5:end]
    end

    config = replace(config, ' ' => "" )

    o1 = String[]
    o2 = String[]
    
    for k ∈ eachindex(nl)
        if occursin(nl[k], config)
            n, ℓ = get(dictAtomicOrbital, nl[k], nothing)
            wl = 2(2ℓ + 1)
            N = wl
            while N > 0
                shl = nl[k] * sup(N)
                if occursin(shl, config)
                    N == wl ? push!(o1, shl) : push!(o2, shl)
                    config = replace(config, shl => "")
                    N = 0
                else
                    N -=1
                end
            end
        end
    end

    iszero(length(config)) || Err && error("Error: $(config) (unknown configuration)") 
    
    return (o1, o2)
        
end

# ------------------------------------------------------------------------------
#                       listConfigurations(Q:Int)
# ------------------------------------------------------------------------------

@doc raw"""
    listConfigurations(Z1::Int, Z2::Int; Q=0)

List of groundstate configurations of atoms with atomic charge `Q`.
#### Example
```
julia> listConfigurations(1:4; Q=0)
(1, " 1s¹", "1s¹")
(2, "1s² ", "[He]")
(3, "1s² 2s¹", "[He]2s¹")
(4, "1s²2s² ", "[Be]")
```
"""
function listConfigurations(Z1::Int, Z2::Int; Q=0)

    a = listAtoms(Z1:Z2, Q; fmt=Object)

    b = unique([(a[i].Z, formatConfig(a[i].config; Err=false), a[i].config) for i ∈ eachindex(a)])

    println.(b)

    return nothing

end
function listConfigurations(itr::UnitRange; Q=0::Int)

    listConfigurations(itr.start, itr.stop; Q)

end

# ------------------------------------------------------------------------------
#                collectConfig(config::String)
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
    config = formatConfig(config)  
    config = replace(config, ' ' => "" )

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