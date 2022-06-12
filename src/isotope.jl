# ======== Isotope(Z, A, radius, mass, I, π, lifetime, mdm, eqm, ra) ===========

"""
    Isotope(symbol Z, A, radius, mass, I, π, lifetime, mdm, eqm, ra)

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
     Z::Int            # atomic number
     N::Int            # neutron number
     A::Int            # atomic mass number (amu)
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

function _specsIsotope(Z::Int, A::Int, isotope)              # Isotope properties

    (symbol, radius, mass, I, π, t½, mdm, eqm, ra) = isotope
    (name, symbol, weight) = get(dictElements, Z, nothing)

    str½ = t½ == 1e100 ? "stable" : "$(t½) years"
    strI = typeof(I) == Rational{Int} ?
          (repr(numerator(I)) * "/" * repr(denominator(I))) : repr(I)

    strRA = ra ≠ nothing ? "$(ra) %" : "trace"
    strmdm = mdm ≠ nothing ? "$(mdm) " : "not available"
    streqm = eqm ≠ nothing ? "$(eqm)" : "not available"

    strIsotope = sup(A) * symbol

    name = (Z,A) == (1,1) ? "hydrogen"  :
           (Z,A) == (1,2) ? "deuterium" :
           (Z,A) == (1,3) ? "tritium" : name

    str = "Isotope created: " * strIsotope * "
    element: " * name * "
    atomic number: Z = $Z
    neutron number: n = $(A-Z)
    atomic mass number: A = $A amu
    rms nuclear charge radius: R = " * repr(radius) * " fm
    atomic mass: mass = " * repr(mass) * " amu
    nuclear spin: I = " * strI * " ħ
    parity of nuclear state: π = $π
    lifetime: " * str½ * "
    nuclear magnetic dipole moment: mdm = " * strmdm * "
    nuclear electric quadrupole moment: eqm = " * streqm * "
    relative abundance: RA = " * strRA

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
* `     .lt`:  lifetime in years (`::Float64`)
#### Examples:
```
isotope = castIsotope(Z=1, A=3, msg=false)
  Isotope(1,-2, 3, 1.7591, 3.016049281, 1//2, 1, 12.33, 2.97896246, 0, nothing)

isotope.ra
  99.9855

castIsotope(Z=1,A=3);
  Isotope created: ³H
      element: tritium
      atomic number: Z = 1
      neutron number: n = 2
      atomic mass number: A =  3 amu
      rms nuclear charge radius: R = 1.7591 fm
      atomic mass: mass = 3.016049281 amu
      nuclear spin: I = 1//2 ħ
      parity of nuclear state: π = 1
      lifetime: 12.33 years
      nuclear magnetic dipole moment: mdm = 2.97896246
      nuclear electric quadrupole moment: eqm = 0
      relative abundance: RA = trace
```
"""
function castIsotope(;Z=1, A=1, msg=true)

    dict = dictIsotopes
    isotope = (Z, A) ∈ keys(dict) ? get(dict, (Z, A), nothing) :
    error("Error: isotope (Z = $Z, A = $A) not present in `dictIsotopes`")

    msg && println(_specsIsotope(Z, A, isotope) )

    symbol, radius, mass, I, π, t½, mdm, eqm, ra = isotope
    symbol = sup(A) * symbol

    return Isotope(symbol, Z, A-Z, A, radius, mass, I, π, t½, mdm, eqm, ra)

end
