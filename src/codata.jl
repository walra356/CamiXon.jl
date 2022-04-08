# ============================== Value(val, unit) ==============================

"""
    Value(val::Real, unit::String)

Object to hold a real numerical value together with a unit specifier.

The fields are:
* ` .val`: numerical value (`::Real`)
* `.unit`: unit specifier (`::String`)
#### Example:
```
f = Value(1,"Hz")
  Value(1, "Hz")

f.val
  1

f.unit
  "Hz"
```
"""
struct Value

    val::Number
    unit::String

end

# ========================== strValue(val) =====================================

"""
    strValue(f::Value)

String expression for a [`Value`](@ref) object in `:compact => true`
representation
#### Example:
```
f = Value(1,"Hz")
strValue(f)
  "1 Hz"
```
"""
function strValue(f::Value)

    strval = repr(f.val, context=:compact => true)
    strunit = " " * f.unit

    return strval * strunit

end

# ============================== NamedValue(val, unit) =========================

"""
    NamedValue(val::Value, name::String, comment::String)

Object to hold a [`Value`](@ref) together with its `symbolic name` and a `short`
description

The fields are:
* `.val`: Value  (`::Value`)
* `.name`: symbolic name (`::String`)
* `.comment`: description (`::String`)

Named Value object
The object `NamedValue` is best created using [`castNamedValue`](@ref).
#### Example:
```
f = Value(1,"Hz")
  Value(1, "Hz", "frequency")

f.name
  "frequency"
```
"""
struct NamedValue

    val::Value
    name::String
    comment::String

end

# ==================== NamedValue(val, unit) ===================================

"""
    castNamedValue(val::Value; name=" ", comment=" ")

Method to create a [`NamedValue`](@ref) object
#### Example
```
v = Value(1.602176634e-19, "C")
nv = castNamedValue(v; name="e")
nv.name * " = " * strValue2(nv.val)
  "e = 1.60218e-19 C"
```
"""
function castNamedValue(val::Value; name=" ", comment=" ")

    return NamedValue(val, name, comment)

end

# ========================= codata =============================================

"""
    Codata

Object to hold the natural constants from CODATA. It is best created with the
function [`castCodata`](@ref)

The fields are:
* `.∆νCs`: Cs hyperfine transition frequency (`::Value`)
* `   .c`: speed of light in vacuum (`::Value`)
* `   .h`: Planck constant (`::Value`)
* `   .ħ`: Planck constant - reduced (`::Value`)
* `   .e`: elementary charge (`::Value`)
* `  .kB`: Boltzmann constant (`::Value`)
* `  .NA`: Avogadro constant (`::Value`)
* ` .Kcd`: Luminous efficacy (`::Value`)
* `  .me`: electron rest mass (`::Value`)
* `  .R∞`: Rydberg constant (`::Value`)
* `  .Ry`: Rydberg frequency (`::Value`)
* `  .Eh`: Hartree a.u. (`::Value`)
* `   .α`: fine-structure constant (`::Value`)
* `  .μ0`: magnetic permitivity of vacuum (`::Value`)
* `  .ε0`: electric permitivity of vacuum (`::Value`)
* `  .KJ`: Josephson constant (`::Value`)
* `  .RK`: Von Klitzing constant (`::Value`)
* `   .R`: Molar gas constant (`::Value`)
* `.matE`: unit conversion matrix (Matrix{Float64})

#### Example:
```
codata = castCodata(2018)
codata.μ0
  Value(1.2566370621250601e-6, "N A⁻²")

codata.μ0.val
  1.2566370621250601e-6
```
"""
struct Codata

    ∆νCs::Value
       c::Value
       h::Value
       ħ::Value
       e::Value
      kB::Value
      NA::Value
     Kcd::Value
      me::Value
      R∞::Value
      Ry::Value
      Eh::Value
       α::Value
      μ0::Value
      ε0::Value
      KJ::Value
      RK::Value
       R::Value
    matE::Matrix{Float64}

end

# ========================= createCodata(year) =================================

"""
    castCodata(year::Int)

Method to create the [`Codata`](@ref) object
#### Example:
```
codata = castCodata(2018)
strValue.([codata.∆νCs,codata.c,codata.h])
 3-element Vector{String}:
  "9192631770 Hz"
  "299792458 m s⁻¹"
  "6.62607e-34 J Hz⁻¹"
```
"""
function castCodata(year::Int)

    year == 2018 || error("Error: codata$(year) not implemented")

    ∆νCs = Value(9192631770, "Hz")
       c = Value(299792458, "m s"*sup(-1))
       h = Value(6.62607015e-34, "J Hz"*sup(-1))
       ħ = Value(h.val/(2π), "J s")
       e = Value(1.602176634e-19, "C")
      kB = Value(1.380649e-23, "J K"*sup(-1))
      NA = Value(6.02214076e23, "mol"*sup(-1))
     Kcd = Value(683, "lm W"*sup(-1))
      me = Value(9.1093837015e-31, "Kg")
      R∞ = Value(10973731.568160, "m"*sup(-1))
      Ry = Value(R∞.val*c.val, "Hz")
      Eh = Value(2Ry.val*h.val, "Hartree a.u.")
       α = Value(sqrt(Eh.val/me.val)/c.val, "")
      μ0 = Value(2α.val*h.val/e.val/e.val/c.val, "N A"*sup(-2))
      ε0 = Value(1/μ0.val/c.val/c.val, "F m"*sup(-1))
      KJ = Value(2e.val/h.val, "Hz V"*sup(-1))
      RK = Value(h.val/e.val/e.val, "Ω")
       R = Value(NA.val*kB.val, "J mol"*sup(-1)*"K"*sup(-1))

    H = 2e-15 * c.val * R∞.val
    J = 1e+15 * h.val
    V = 1e+15 * h.val / e.val

    M =[1 1e3 1e6 1e9 1e12 1e15 1e18 1e21 1e24 1e+21*H 0.5e+21*H 1e21/J 1e21/V;
        1 1 1e3 1e6 1e9 1e12 1e15 1e18 1e21 1e+18*H 0.5e+18*H 1e18/J 1e18/V;
        1 1 1 1e3 1e6 1e9 1e12 1e15 1e18 1e+15*H 0.5e+15*H 1e15/J 1e15/V;
        1 1 1 1 1e3 1e6 1e9 1e12 1e15 1e+12*H 0.5e+12*H 1e12/J 1e12/V;
        1 1 1 1 1 1e3 1e6 1e9 1e12 1e+9*H 0.5e+9*H 1e9/J 1e9/V;
        1 1 1 1 1 1 1e3 1e6 1e9 1e+6*H 0.5e+6*H 1e6/J 1e6/V;
        1 1 1 1 1 1 1 1e3 1e6 1e3*H 0.5e+3*H 1e3/J 1e3/V;
        1 1 1 1 1 1 1 1 1e3 1*H 0.5*H 1/J 1/V;
        1 1 1 1 1 1 1 1 1 1e-3*H 0.5e-3*H 1e-3/J 1e-3/V;
        1 1 1 1 1 1 1 1 1 1 0.5 2.2937122783963e17 1/27.211386245988;
        1 1 1 1 1 1 1 1 1 1 1 2*2.2937122783963e17 2/27.211386245988;
        1 1 1 1 1 1 1 1 1 1 1 1 1/6.241509074e18;
        1 1 1 1 1 1 1 1 1 1 1 1 1]

    N = Int(sqrt(length(M)))

    for i=1:N
        for j=i:N
            M[j,i] = 1/M[i,j]
        end
    end

    o = Codata(∆νCs,c,h,ħ,e,kB,NA,Kcd,me,R∞,Ry,Eh,α,μ0,ε0,KJ,RK,R,M)

    return o

end

# ========================= listCodata(codata::Codata) =========================

"""
    listCodata(codata::Codata)

Method to list the fields of [`Codata`](@ref) by their symbolic name
#### Example:
```
codata = castCodata(2018)
listCodata(codata::Codata)
  ∆νCs = 9192631770 Hz
  c = 299792458 m s⁻¹
  h = 6.62607e-34 J Hz⁻¹
  ħ = 1.05457e-34 J s
  e = 1.60218e-19 C
  kB = 1.38065e-23 J K⁻¹
  NA = 6.02214e23 mol⁻¹
  Kcd = 683 lm W⁻¹
  mₑ = 9.10938e-31 Kg
  R∞ = 1.09737e7 m⁻¹
  Ry = 3.28984e15 Hz
  Eₕ = 4.35974e-18 Hartree a.u.
  α = 0.00729735
  μ₀ = 1.25664e-6 N A⁻²
  ε₀ = 8.85419e-12 F m⁻¹
  KJ = 4.83598e14 Hz V⁻¹
  RK = 25812.8 Ω
  R = 8.31446 J mol⁻¹K⁻¹
```
"""
function listCodata(codata::Codata)

    ∆νCs = NamedValue(codata.∆νCs, "∆νCs",
                       sup(133)*"Cs hyperfine transition frequency")
       c = NamedValue(codata.c, "c", "speed of light in vacuum")
       h = NamedValue(codata.h, "h", "Planck constant")
       ħ = NamedValue(codata.ħ, "ħ", "Planck constant (reduced)")
       e = NamedValue(codata.e, "e", "elementary charge")
      kB = NamedValue(codata.kB, "kB", "Boltzmann constant")
      NA = NamedValue(codata.NA, "NA", "Avogadro constant")
     Kcd = NamedValue(codata.Kcd, "Kcd", "Luminous efficacy")
      me = NamedValue(codata.me, "m"*sub("e"), "electron rest mass")
      R∞ = NamedValue(codata.R∞, "R∞", "Rydberg constant")
      Ry = NamedValue(codata.Ry, "Ry", "Rydberg frequency")
      Eh = NamedValue(codata.Eh, "E"*sub("h"), "Hartree atomic unit")
       α = NamedValue(codata.α, "α", "fine-structure constant")
      μ0 = NamedValue(codata.μ0, "μ"*sub(0), "magnetic permitivity of vacuum")
      ε0 = NamedValue(codata.ε0, "ε"*sub(0), "electric permitivity of vacuum")
      KJ = NamedValue(codata.KJ, "KJ", "Josephson constant")
      RK = NamedValue(codata.RK, "RK", "Von Klitzing constant")
       R = NamedValue(codata.R, "R", "Molar gas constant")

    U =  [∆νCs, c, h, ħ, e, kB, NA, Kcd, me, R∞, Ry, Eh, α, μ0, ε0, KJ, RK, R]

    for i ∈ eachindex(U)
        println(U[i].name * " = " * strValue(U[i].val))
    end

    return

end

# =============== convertUnit(val; unitIn="kHz", unitOut="xHz") ===============

"""
    convertUnit(val, codata::Codata; unitIn="Hartree", unitOut="xHz")

Unit conversion between μHz, ..., EHz, Hartree, Rydberg, Joule, and eV

default input: Hartree

default output: xHz ∈ {μHz, mHz, Hz, kHz, MHz, GHz, THz, PHz, EHz}
#### Example:
```
codata = createCodata(2018)
convertUnit(1, codata; unitIn="Hz", unitOut="Joule")
  6.62607015e-34

convertUnit(1, codata; unitIn="Hartree", unitOut="Hz")
  Value(6.57968392050182e15, "Hz")

f = convertUnit(1, codata) # default input (Hartree) and output (xHz)
strf = strValue(f)
  "6.57968 PHz"
```
"""
function convertUnit(val, codata::Codata; unitIn="Hartree", unitOut="xHz")
# ==============================================================================
#  convertUnit(val; unitIn="kHz", unitOut="MHz")
# ==============================================================================
    U = ["μHz","mHz","Hz","kHz","MHz","GHz","THz","PHz","EHz"]
        push!(U,"Hartree")
        push!(U,"Rydberg")
        push!(U,"Joule")
        push!(U,"eV")

    unitIn ∈ U || error("Error: unitIn = $(unitIn) unknown unit type")

    M = codata.matE

    v = zeros(Float64,length(U))
    i = findfirst(x -> x==unitIn, U)
    v[i] = val
    w = M * v

    if unitOut == "xHz"
        i = findfirst(x -> x=="Hz", U)
        unitOut = 1e18 ≤  w[i] ? "EHz" :
                  1e15 ≤  w[i] < 1e18 ? "PHz" :
                  1e12 ≤  w[i] < 1e15 ? "THz" :
                  1e9  ≤  w[i] < 1e12 ? "GHz" :
                  1e6  ≤  w[i] < 1e9  ? "MHz" :
                  1e3  ≤  w[i] < 1e6  ? "kHz" :
                  1e0  ≤  w[i] < 1e3  ? "Hz"  :
                  1e-3 ≤  w[i] < 1e0  ? "mHz" : "μHz"
    else
        unitOut ∈ U || error("Error: unitOut = $(unitOut) unknown unit type")
    end

    i = findfirst(x -> x==unitOut, U)

    return Value(w[i],unitOut)

end

# ============ calibrationReport(E, Ecal; unitIn="Hartree") ====================

"""
    calibrationReport(E, Ecal, codata::Codata; unitIn="Hartree")

Comparison of energy E with calibration value Ecal

default input: Hartree
#### Example:
```
codata = castCodata(2018)
calibrationReport1(1.1, 1.0, codata; unitIn="Hartree")
  calibration report (Float64):
  Ecal = 1.0 Hartree
  E = 1.1 Hartree
  absolute accuracy: ΔE = 0.1 Hartree (657.968 THz)
  relative accuracy: ΔE/E = 0.0909091
```
"""
function calibrationReport(E, Ecal, codata::Codata; unitIn="Hartree")

    T = typeof(E)

    V = BigFloat

    E = myconvert(V, E)
    Ecal = myconvert(V, Ecal)

    ΔE = abs(E-Ecal)
    ΔErel = ΔE/E

    Δf = convertUnit(ΔE, codata)
    strΔf = " (" * strValue(Δf) * ")"
    strΔE = repr(ΔE, context=:compact => true)
    strΔErel = repr(ΔErel, context=:compact => true)

    msg = "calibration report ($T):\n"
    msg *= @sprintf "Ecal = %.17g %s \n" Ecal unitIn
    msg *= @sprintf "E = %.17g %s \n" E unitIn
    msg *= ΔE ≠ 0 ? "absolute accuracy: ΔE = " * strΔE * " " * unitIn * strΔf * "\n" :
                    "absolute accuracy: ΔE = 0 (exact within $T precision)\n"
    msg *= ΔE ≠ 0 ? "relative accuracy: ΔE/E = " * strΔErel * "\n"                   :
                    "relative accuracy: ΔE/E = 0 (exact within $T precision)\n"

    return println(msg)

end

# ==================== myconvert(T::Type, val::V) ==============================

@doc raw"""
    myconvert(T::Type, val::V) where V <: Number

Conversion including `BigFloat`, `Rational{BigFloat}`, BigInt, and
`Rational{BigInt}`
#### Examples:
```
convert(BigInt,1//3)
  InexactError: BigInt(1//3)

myconvert(BigInt, 1//3)
  0.3333333333333333333333333333333333333333333333333333333333333333333333333348
```
"""
function myconvert(T::Type, val::V) where V <: Number      #### moet verplaatst?
# ==============================================================================
# myconvert(T::Type, val::V) # generalization of convert to include BigFloat
# ==============================================================================

    if T == BigFloat
        o = V <: Rational ?
            T(string(numerator(val)))/T(string(denominator(val))) :
            T(string(val))
    elseif T == BigInt
        o = V <: Rational ?
            T(numerator(val))/T(denominator(val)) :
            T(val)
    else
        o = convert(T,val)
    end

    return o

end
