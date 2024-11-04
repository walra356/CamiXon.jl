# SPDX-License-Identifier: MIT

# author: Jook Walraven - 14-2-2023

# ==============================================================================
#                               thermal-properties.jl
# ==============================================================================

@doc raw"""
    melting_point(atomicnumber::Int)
    melting_point(element::String)

Melting point of a given `element` at standard pressure (1 atm).

#### Example:
```
julia> melting_point("Li")
453.65
```
"""
function melting_point(atomicnumber::Int)

    return dictMeltingPoints[atomicnumber]

end
function melting_point(element::String)

    A = dictAtomicNumbers[element]

    return melting_point(A)

end

@doc raw"""
    svp(atomicnumber::Int, temp::Real)
    svp(element::String, temp::Real)

Saturated vapor pressure, *p* (in Pa), of a given `element` for a given 
temperature *T* (in K),

```math
\mathrm{log_{e}}p=A+B/T+C\mathrm{log_{10}}T+D\cdot T/1000,
```
where A,B,C,D, are the Antoine coefficients collected in 
[`CamiXon.dictAntoineCoefficients`](@ref).

#### Examples:
```
julia> svp("Li", 623.0)
0.0015230367024569058
```
"""
function svp(atomicnumber::Int, temp::Real)

    T = float(temp)
    i = atomicnumber
    d = CamiXon.dictAntoineCoefficients

    groupA = [3, 4, 11, 13, 19, 21, 22, 23, 26, 27, 28, 29, 30, 31, 37, 39, 40, 45, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 64, 65, 68, 71, 78, 79, 80, 81, 82, 90, 91, 92, 93, 94, 96]
    groupB = [12, 20, 24, 25, 38, 41, 42, 44, 62, 63, 66, 67, 69, 70, 72, 73, 74, 75, 76, 77, 95]
    groupC = [1, 2, 5, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36, 43, 51, 52, 53, 54, 61, 83, 84, 85, 86, 87, 88, 89, 97, 98, 99, 100, 101, 102]
    groupD = (6)

    (Tmin, mp, Tmax) = d[i][3]

    if i ∈ groupA
        (A, B, C, D) = (T ≤ mp) ? d[i][1] : d[i][2]
    elseif i ∈ groupB
        (A, B, C, D) = (T ≤ mp) ? d[i][1] :
            println("Antoine coefficients for T > m.p. = $(mp) not found")
    else
        println("Antoine coefficients not found")
    end

    lp = A + B / T + C * log10(T) + D * T / 1000.0

    p = exp(lp) # saturated vapor pressure

    return p

end
function svp(element::String, temp::Real)

    A = dictAtomicNumbers[element]
    p = svp(A, temp)

    return p

end

@doc raw"""
    latent_heat_vaporization(atomicnumber::Int, temp::Real)
    latent_heat_vaporization(element::String, temp:Real)

Latent heat of vaporization, *L*(*T*) (in Joule/K), of a given `element` 
for temperature *T* (in K), 

```math
L(T) = -(B +C\cdot T \mathrm{log_{10}}T+D\cdot T^2/1000),
```

where A,B,C,D, are the Antoine coefficients collected in 
[`CamiXon.dictAntoineCoefficients`](@ref).

#### Example:
```
julia> latent_heat_vaporization("Li", 623.0)
-18473.64020109123
```
"""
function latent_heat_vaporization(atomicnumber::Int, temp::Real)

    T = float(temp)
    i = atomicnumber
    d = CamiXon.dictAntoineCoefficients

    groupA = [3, 4, 11, 13, 19, 21, 22, 23, 26, 27, 28, 29, 30, 31, 37, 39, 40, 45, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 64, 65, 68, 71, 78, 79, 80, 81, 82, 90, 91, 92, 93, 94, 96]
    groupB = [12, 20, 24, 25, 38, 41, 42, 44, 62, 63, 66, 67, 69, 70, 72, 73, 74, 75, 76, 77, 95]
    groupC = [1, 2, 5, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36, 43, 51, 52, 53, 54, 61, 83, 84, 85, 86, 87, 88, 89, 97, 98, 99, 100, 101, 102]
    groupD = (6)

    (Tmin, mp, Tmax) = d[i][3]

    if i ∈ groupA
        (A, B, C, D) = (T ≤ mp) ? d[i][1] : d[i][2]
    elseif i ∈ groupB
        (A, B, C, D) = (T ≤ mp) ? d[i][1] :
                println("Antoine coefficients for T > m.p. = $(mp) not found")
    else
        println("Antoine coefficients not found")
    end

    L = B + C * T * log10(T) + D * T * T / 1000.0 # latent heat 0f vaporization

    return -L

end
function latent_heat_vaporization(element::String, temp::Real)

    atomicnumber = dictAtomicNumbers[element]

    return latent_heat_vaporization(atomicnumber, temp)

end