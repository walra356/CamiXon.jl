# SPDX-License-Identifier: MIT

# author: Jook Walraven - 14-2-2023

# ==============================================================================
#                               isotope.jl
# ==============================================================================

# ======== Isotope(Z, A, radius, mass, I, π, lifetime, mdm, eqm, ra) ===========

"""
    Isotope(symbol, name, Z, A, N, R, M, I, π, T½, mdm, eqm, ra)

Type with fields:
* `     .symbol`: symbol (`::String`)
* `     .name`: name (`::String`)
* `     .Z`:  atomic number (`::Int`)
* `     .A`:  atomic mass number in amu (`::Int`)
* `     .N`:  neutron number (`::Int`)
* `     .R`:  rms charge radius in Fermi (`::Float64`)
* `     .M`:  *atomic* mass in amu (`::Float64`)
* `     .I`:  nuclear spin in units of ħ  (`::Rational{Int}`)
* `     .π`:  parity of nuclear state (`::Int`)
* `     .T½`:  lifetime in years (`::Float64`)
* `     .mdm`: nuclear magnetic dipole moment (`::Float64`)
* `     .eqm`: nuclear electric quadrupole moment (`::Float64`)
* `     .ra`:  relative abundance in % (`::Float64`)

The type `Isotope` is best created with the function [`castIsotope`](@ref).
"""
struct Isotope                     # Isotopic properties
     symbol::String                # isotope symbol
     name::String                  # name isotope
     Z::Int                        # atomic number
     A::Int                        # atomic mass number (amu)
     N::Int                        # neutron number
     R::Float64                    # rms charge radius (Fermi)
     M::Float64                    # nuclear mass (amu)
     I::Union{Rational{Int}, Int}  # nuclear spin in units of ħ
     π::Int                        # parity of nuclear state
     T½::Float64                   # lifetime (years)
     mdm::Float64                  # nuclear magnetic dipole moment
     eqm::Union{Float64, Nothing}  # nuclear electric quadrupole moment
     ra::Union{Float64, Nothing}   # relative abundance (%)
end

# ..............................................................................
function _stdIsotope(Z::Int, A::Int)

    dict = dictIsotopes

    o = (Z, A) ∈ keys(dict) ? castIsotope(;Z, A, msg=false) : return nothing

    return o

end
# ..............................................................................
function _strIsotope(Z::Int, A::Int)

    dict = dictIsotopes
    isotope = (Z, A) ∈ keys(dict) ? castIsotope(;Z, A, msg=false) : return nothing

    strπ = isotope.π == 1 ? "⁺" : "⁻"
    strRA = isnothing(isotope.ra) ? "trace" : repr(isotope.ra) * "%"
    strT½ = isotope.T½ == 1e100 ? "stable" : "radioactive"

    str = isotope.symbol
    str *= ", " * isotope.name
    str *= ", Z=" * repr(isotope.Z)
    str *= ", A=" * repr(isotope.A)
    str *= ", N=" * repr(isotope.N)
    str *= ", R=" * repr(isotope.R)
    str *= ", M=" * repr(isotope.M)
    str *= ", I=" * strRational(isotope.I) * strπ
    str *= ", μI=" * repr(isotope.mdm)
    str *= ", Q=" * repr(isotope.eqm)
    str *= ", RA=" * strRA
    str *= ", (" * strT½ * ")"

    return str

end
# ..............................................................................
function _texIsotope(Z::Int, A::Int; indent=false)              # Isotope properties

    dict = dictIsotopes
    isotope = (Z, A) ∈ keys(dict) ? castIsotope(;Z, A, msg=false) : return nothing
    (name, symbol, weight) = get(dictElements, Z, nothing)

    strπ = isotope.π == 1 ? "\$^+\$" : "\$^-\$"
    name = isotope.name
    strRA = isnothing(isotope.ra) ? "trace" : repr(isotope.ra)
    strT½ = isotope.T½ == 1e100 ? "\\," : "*\$\\!\\!\$"
    symbol = name=="deuterium" ? "D" : name=="tritium" ? "T" : symbol

    str = indent ? "" : repr(isotope.Z)
    str *= " & " * (!indent ? name : name=="deuterium" ? name : name=="tritium" ? name : "")
    str *= " & " * "\$^{$A}\$" * symbol
    str *= " & " * repr(isotope.A) * strT½
    str *= " & " * repr(isotope.N)
    str *= " & " * repr(isotope.R)
    str *= " & " * repr(isotope.M)
    str *= " & " * strRational(isotope.I) * strπ
    str *= " & " * repr(isotope.mdm)
    str *= " & " * repr(isotope.eqm)
    str *= " & " * strRA
    str *= " \\\\\n"

    return str

end
# ..............................................................................
function _infoIsotope(Z::Int, A::Int)

    dict = dictIsotopes
    isotope = (Z, A) ∈ keys(dict) ? castIsotope(;Z, A, msg=false) : return nothing

    strπ = isotope.π == 1 ? "even" : "odd"
    strRA = isnothing(isotope.ra) ? "trace" : repr(isotope.ra) * "%"
    strT½ = isotope.T½ == 1e100 ? "stable" : repr(isotope.T½) * " years"

    str = "Isotope: " * isotope.name * "-" * repr(isotope.A)
    str *= "\n  symbol: " * isotope.symbol
    str *= "\n  element: " * isotope.name
    str *= "\n  atomic number: Z = " * repr(isotope.Z)
    str *= "\n  atomic mass number: A = " * repr(isotope.A)
    str *= "\n  neutron number: N = " * repr(isotope.N)
    str *= "\n  rms nuclear charge radius: R = " * repr(isotope.R) * " fm"
    str *= "\n  atomic mass: M = " * repr(isotope.M) * " amu"
    str *= "\n  nuclear spin: I = " * strRational(isotope.I) * " ħ"
    str *= "\n  parity of nuclear state: π = " * strπ
    str *= "\n  nuclear magnetic dipole moment: μI = " * repr(isotope.mdm) * " μN"
    str *= "\n  nuclear electric quadrupole moment: Q = " * repr(isotope.eqm) * " barn"
    str *= "\n  relative abundance: RA = " * strRA
    str *= "\n  lifetime: " * strT½

    return println(str)

end
# ..............................................................................
"""
    listIsotope(Z::Int, A::Int; fmt=Object)

Properties of isotopes with atomic number `Z` and atomic mass number `A`.

Output options: `fmt` =  `Object` (default), `String`, `Latex`, `Info`.
#### Example:
```
listIsotope(1,3; fmt=Info)
  Isotope: tritium-3
    symbol: ³T
    element: tritium
    atomic number: Z = 1
    atomic mass number: A = 3
    neutron number: N = 2
    rms nuclear charge radius: R = 1.7591 fm
    atomic mass: M = 3.016049281 amu
    nuclear spin: I = 1/2 ħ
    parity of nuclear state: π = even
    nuclear magnetic dipole moment: μI = 2.97896246μN
    nuclear electric quadrupole moment: Q = 0.0barn
    relative abundance: RA = trace
    lifetime: 12.33 years
```
"""
function listIsotope(Z::Int, A::Int; fmt=Object)

    fmt === Object && return _stdIsotope(Z, A)
    fmt === String && return _strIsotope(Z, A)
    fmt === Latex && return _texIsotope(Z, A)
    fmt === Info && return _infoIsotope(Z, A)

    return error("Error: unknown output type")

end
# ..............................................................................
"""
    listIsotopes(Z1::Int, Z2::Int; fmt=Object)

All isotopes with atomic number from `Z1` to `Z2`.

Output options: `Object` (default), `String`, `Latex`, `Info`.
#### Example:
```
listIsotopes(1,3) == listIsotopes(1:3)
 true

listIsotopes(1:1; fmt=Info)
 3-element Vector{Any}:
  Isotope("¹H", "hydrogen", 1, 1, 0, 0.8783, 1.007825032, 1//2, 1, 1.0e100, 2.792847351, 0.0, 99.9855)
  Isotope("²D", "deuterium", 1, 2, 1, 2.1421, 2.014101778, 1, 1, 1.0e100, 0.857438231, 0.0028578, 0.0145)
  Isotope("³T", "tritium", 1, 3, 2, 1.7591, 3.016049281, 1//2, 1, 12.33, 2.97896246, 0.0, nothing)
```
"""
function listIsotopes(Z1::Int, Z2::Int; fmt=Object)

    o = []

    for Z=Z1:Z2
        for A=1:3Z
            next = listIsotope(Z, A; fmt)
            isnothing(next) ? false : Base.push!(o, next)
        end
    end

    return o

end
function listIsotopes(itrZ; fmt=Object)

    return listIsotopes(itrZ.start, itrZ.stop; fmt)

end

# ================= castIsotope(;Z=1, A=1, msg=true) ===========================

#...............................................................................
function _filter_isotopes(elt::String, A::Int)

    elt == "D" && return  A == 2 ? 2 :
                          error("Error: A = 2 expected for isotope $(elt) ")
    elt == "T" && return  A == 3 ? 3 :
                          error("Error: A = 3 expected for isotope $(elt) ")

    return A

end
#...............................................................................
"""
    castIsotope(;Z=1, A=1, msg=false)
    castIsotope(elt::String; A=1, msg=false)

Create Isotope with fields
* `     .symbol`: symbol (`::String`)
* `     .name`: symbol (`::String`)
* `     .Z`:  atomic number (`::Int`)
* `     .A`:  atomic mass number in amu (`::Int`)
* `     .N`:  neutron number (`::Int`)
* `     .R`:  rms charge radius in Fermi (`::Float64`)
* `     .M`:  atomic mass in amu (`::Float64`)
* `     .I`:  nuclear spin in units of ħ (`::Rational{Int}`)
* `     .π`:  parity of nuclear state (`::Int`)
* `     .ra`:  relative abundance in % (`::Float64`)
* `     .mdm`: nuclear magnetic dipole moment (`::Float64`)
* `     .eqm`: nuclear electric quadrupole moment (`::Float64`)
* `     .T½`:  lifetime in years (`::Float64`)
#### Examples:
```
castIsotope("Rb"; A=87) == castIsotope(Z=37, A=87)
  true

isotope = castIsotope(Z=1, A=3)
  Isotope("³T", "tritium", 1, 3, 2, 1.7591, 3.016049281, 1//2, 1, 12.33, 2.97896246, 0, nothing)

isotope.T½
  12.33

castIsotope(Z=1, A=3, msg=true);
  Isotope created: tritium-3
      symbol: ³T
      element: tritium
      atomic number: Z = 1
      atomic mass number: A = 3
      neutron number: N = 2
      rms nuclear charge radius: R = 1.7591 fm
      atomic mass: M = 3.016049281 amu
      nuclear spin: I = 1/2 ħ
      parity of nuclear state: π = ⁺
      nuclear magnetic dipole moment: μI = 2.97896246μN
      nuclear electric quadrupole moment: Q = 0.0barn
      relative abundance: RA = trace
      lifetime: 12.33 years
```
"""
function castIsotope(;Z=1, A=1, msg=false)

    dict = dictIsotopes
    isotope = (Z, A) ∈ keys(dict) ? get(dict, (Z, A), nothing) :
    error("Error: isotope (Z = $Z, A = $A) not present in `dictIsotopes`")

    (symbol, name, Z, A, N, R, M, I, π, T½, mdm, eqm, ra) = isotope

    o = Isotope(symbol, name, Z, A, N, R, M, I, π, T½, mdm, eqm, ra)

    msg && println("Isotope created: " * listIsotope(Z, A ; fmt=String))

    return o

end
function castIsotope(elt::String; A=1, msg=false)

    dict = dictAtomicNumbers
    Z = (elt) ∈ keys(dict) ? get(dict, elt, nothing) :
        return error("Error: element $(elt) - not found in `dictAtomNumbers`")

    A = _filter_isotopes(elt, A)

    dict = dictIsotopes
    isotope = (Z, A) ∈ keys(dict) ? get(dict, (Z, A), nothing) :
    error("Error: isotope (Z = $Z, A = $A) not present in `dictIsotopes`")

    (symbol, name, Z, A, N, R, M, I, π, T½, mdm, eqm, ra) = isotope

    o = Isotope(symbol, name, Z, A, N, R, M, I, π, T½, mdm, eqm, ra)

    msg && println("Isotope created: " * listIsotope(Z, A ; fmt=String))

    return o

end

# ============================== end ===========================================
