# ======== Isotope(Z, A, radius, mass, I, π, lifetime, mdm, eqm, ra) ===========

"""
    Isotope(symbol, name, Z, A, N, R, M, I, π, t½, mdm, eqm, ra)

Type with fields:
* `     .symbol`: symbol (`::String`)
* `     .Z`:  atomic number (`::Int`)
* `     .N`:  neutron number (`::Int`)
* `     .A`:  atomic mass number in amu (`::Int`)
* `     .R`:  rms charge radius in Fermi (`::Float64`)
* `     .M`:  *atomic* mass in amu (`::Float64`)
* `     .I`:  nuclear spin in units of ħ  (`::Rational{Int}`)
* `     .π`:  parity of nuclear state (`::Int`)
* `     .t½`:  lifetime in years (`::Float64`)
* `     .mdm`: nuclear magnetic dipole moment (`::Float64`)
* `     .eqm`: nuclear electric quadrupole moment (`::Float64`)
* `     .ra`:  relative abundance in % (`::Float64`)

The type `Isotope` is best created with the function [`castIsotope`](@ref).
"""
struct Isotope         # Isotopic properties
     symbol::String    # isotope symbol
     name::String      # name isotope
     Z::Int            # atomic number
     A::Int            # atomic mass number (amu)
     N::Int            # neutron number
     R::Float64        # rms charge radius (Fermi)
     M::Float64        # nuclear mass (amu)
     I::Rational{Int}  # nuclear spin in units of ħ
     π::Int            # parity of nuclear state
     t½::Float64       # lifetime (years)
     mdm::Float64      # nuclear magnetic dipole moment
     eqm::Union{Float64, Nothing}       # nuclear electric quadrupole moment
     ra::Union{Float64, Nothing}       # relative abundance (%)
end

# ====================== _specsIsotope(Z, A) ===================================

function _specsIsotope(isotope)

    dict = dictIsotopes
    isotope = (Z, A) ∈ keys(dict) ? castIsotope(;Z, A, msg=false) : return nothing

    strπ = isotope.π == 1 ? "⁺" : "⁻"
    strRA = isotope.ra == nothing ? "trace" : repr(isotope.ra) * "%"
    strt½ = isotope.t½ == 1e100 ? "stable" : "radioactive"

    str = isotope.symbol * "
    element: " * name * "
    atomic number: Z = " * repr(isotope.Z) * "
    atomic mass number: A = " * repr(isotope.A) * "
    neutron number: N = " * repr(isotope.N) * "
    rms nuclear charge radius: R = " * repr(isotope.R) * " fm
    atomic mass: M = " * repr(isotope.M) * " amu
    nuclear spin: I = " * strRational(isotope.I) * " ħ
    parity of nuclear state: π = " * strπ * "
    nuclear magnetic dipole moment: μI = " * repr(isotope.mdm) * "
    nuclear electric quadrupole moment: Q = " * repr(isotope.eqm) * "
    relative abundance: RA = " * strRA * " %
     (" * strt½ * ")"

    return str

end

# ================= castIsotope(;Z=1, A=1, msg=true) ===========================

"""
    castIsotope(;Z=1, A=1, msg=true)

Create Isotope with fields
* `     .symbol`: symbol (`::String`)
* `     .Z`:  atomic number (`::Int`)
* `     .N`:  neutron number (`::Int`)
* `     .A`:  atomic mass number in amu (`::Int`)
* `     .R`:  rms charge radius in Fermi (`::Float64`)
* `     .M`:  atomic mass in amu (`::Float64`)
* `     .I`:  nuclear spin in units of ħ (`::Rational{Int}`)
* `     .π`:  parity of nuclear state (`::Int`)
* `     .ra`:  relative abundance in % (`::Float64`)
* `     .mdm`: nuclear magnetic dipole moment (`::Float64`)
* `     .eqm`: nuclear electric quadrupole moment (`::Float64`)
* `     .t½`:  lifetime in years (`::Float64`)
#### Examples:
```
isotope = castIsotope(Z=1, A=3, msg=false)
  Isotope("³T", tritium, 1, 3, 2, 1.7591, 3.016049281, 1//2, 1, 12.33, 2.97896246, 0, nothing)

isotope.ra
  99.9855

castIsotope(Z=1,A=3);
  Isotope created: ³T
      element: tritium
      atomic number: Z = 1
      atomic mass number: A =  3
      neutron number: N = 2
      rms nuclear charge radius: R = 1.7591 fm
      atomic mass: M = 3.016049281 amu
      nuclear spin: I = 1/2 ħ
      parity of nuclear state: π = 1
      lifetime: 12.33 years
      nuclear magnetic dipole moment: μI = 2.97896246
      nuclear electric quadrupole moment: Q = 0
      relative abundance: RA = trace
```
"""
function castIsotope(;Z=1, A=1, msg=true)

    dict = dictIsotopes
    isotope = (Z, A) ∈ keys(dict) ? get(dict, (Z, A), nothing) :
    error("Error: isotope (Z = $Z, A = $A) not present in `dictIsotopes`")

    msg && println("Isotope created: " * _specsIsotope(isotope) )

    return Isotope(isotope)

end
