# ============================== Value(val, unit) ==============================

"""
    Value(val::Real, unit::String)

Type with fields:
* ` .val`: numerical value
* `.unit`: unit specifier

Value object
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

String expression for Value object in `:compact => true` representation
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

Type with fields:
* `.val`: Value
* `.name`: symbolic name
* `.comment`: description

Named Value object
The object `NamedValue` is best created by the function `createNamedValue`.
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
    createNamedValue(val::Value; name=" ", comment=" ")

Method to create NamedValue object
#### Example
```
v = Value(1.602176634e-19, "C")
nv = createNamedValue(v; name="e")
nv.name * " = " * strValue2(nv.val)
 "e = 1.60218e-19 C"
```
"""
function createNamedValue(val::Value; name=" ", comment=" ")

    return NamedValue(val, name, comment)

end

# ========================= codata =============================================

"""
    Codata

* `.∆νCs`: Cs hyperfine transition frequency
* `   .c`: speed of light in vacuum
* `   .h`: Planck constant
* `   .ħ`: Planck constant (reduced)
* `   .e`: elementary charge
* `  .kB`: Boltzmann constant
* `  .NA`: Avogadro constant
* ` .Kcd`: Luminous efficacy
* `  .me`: electron rest mass
* `  .R∞`: Rydberg constant
* `  .Ry`: Rydberg frequency
* `  .Eh`: Hartree a.u.
* `   .α`: fine-structure constant
* `  .μ0`: magnetic permitivity of vacuum
* `  .ε0`: electric permitivity of vacuum
* `  .KJ`: Josephson constant
* `  .RK`: Von Klitzing constant
* `   .R`: Molar gas constant
* `.matE`: unit conversion matrix

Codata object
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
    createCodata(year::Int)

Create codata object
#### Example:
```
codata = createCodata(2018)
strValue.([codata.∆νCs,codata.c,codata.h])
3-element Vector{String}:
 "9192631770 Hz"
 "299792458 m s⁻¹"
 "6.62607e-34 J Hz⁻¹"
```
"""
function createCodata(year::Int)

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

# ========================= listCodata(codata::Codata) =================================

"""
    listCodata(codata::Codata)

List codata values by name
#### Example:
```
# to be done
```
"""
function listCodata(codata::Codata)

    ∆νCs = NamedValue(codata.∆νCs, "∆νCs", sup(133)*"Cs hyperfine transition frequency")
       c = NamedValue(codata.c, "c", "speed of light in vacuum")
       h = NamedValue(codata.h, "h", "Planck constant")
       ħ = NamedValue(codata.ħ, "ħ", "Planck constant (reduced)")
       e = NamedValue(codata.e, "e", "elementary charge")
      kB = NamedValue(codata.kB, "kB", "Boltzmann constant")
      NA = NamedValue(codata.NA, "NA", "Avogadro constant")
     Kcd = NamedValue(codata.Kcd, "Kcd", "Luminous efficacy", "Kcd")
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
        println(U[i].name * " = " * strVal(U[i].val))
    end

    return

end

# =============== convertUnits(val; unitIn="kHz", unitOut="xHz") ===============

"""
    convertUnits(val, codata::Codata; unitIn="Hartree", unitOut="xHz")

Unit conversion between μHz, ..., EHz, Hartree, Rydberg, Joule, and eV

default input: Hartree

default output: xHz ∈ {μHz, mHz, Hz, kHz, MHz, GHz, THz, PHz, EHz}
#### Example:
```
codata = createCodata(2018)
convertUnits(1, codata; unitIn="Hz", unitOut="Joule")
 6.62607015e-34

convertUnits(1, codata; unitIn="Hartree", unitOut="Hz")
 Value(6.57968392050182e15, "Hz")

f = convertUnits(1, codata) # default input (Hartree) and output (xHz)
strf = strValue(f)
 "6.57968 PHz"
```
"""
function convertUnits(val, codata::Codata; unitIn="Hartree", unitOut="xHz")
# ==============================================================================
#  convertUnits(val; unitIn="kHz", unitOut="MHz")
# ==============================================================================
    U = ["μHz","mHz","Hz","kHz","MHz","GHz","THz","PHz","EHz","Hartree","Rydberg","Joule","eV"]

    unitIn ∈ U || error("Error: unitIn ∉ {μHz, ..., EHz, Hartree, Rydberg, Joule, eV}")

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
        unitOut ∈ U || error("Error: unitOut ∉ {μHz, ..., EHz, Hartree, Rydberg, Joule, eV}")
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
calibrationReport(1.1, 1.0; unitIn="Hartree")
 calibration report:
 absolute accuracy: ΔE = 0.1 Hartree (657.968 THz)
 relative accuracy: ΔE/E = 0.0909091
 Ecal = 1.0 Hartree
 E = 1.100000000000000000000000000000000000000000000000000000000000000000000000000003 Hartree
 input number type: Float64
```
"""
function calibrationReport(E, Ecal, codata::Codata; unitIn="Hartree")

    T = typeof(E)

    V = BigFloat

    E = myconvert(V, E)
    Ecal = myconvert(V, Ecal)

    ΔE = abs(E-Ecal)
    ΔErel = ΔE/E

    Δf = convertUnits(ΔE, codata)
    strΔf = strValue(Δf)
    strΔE = repr(ΔE, context=:compact => true)
    strΔErel = repr(ΔErel, context=:compact => true)

    msg = "calibration report:\n"
    msg *= "absolute accuracy: ΔE = " * strΔE * " " * unitIn * " (" * strΔf * ")\n"
    msg *= "relative accuracy: ΔE/E = " * strΔErel * "\n"
    msg *= "Ecal = $(Ecal) " * unitIn * "\n"
    msg *= "E = $E " * unitIn * "\n"
    msg *= "input number type: $(T)"

    return println(msg)

end

# ==================== myconvert(T::Type, val::V) ==============================

@doc raw"""
    myconvert(T::Type, val::V) where V <: Number

Conversion including BigFloat and BigInt
#### Examples:
```
convert(BigInt,1//3)
 InexactError: BigInt(1//3)

myconvert(BigInt, 1//3)
 0.3333333333333333333333333333333333333333333333333333333333333333333333333333348
```
"""
function myconvert(T::Type, val::V) where V <: Number        #### moet verplaatst?
# ================================================================================
# myconvert(T::Type, val::V) # generalization of convert to include BigFloat
# ================================================================================

    if T == BigFloat
        o = V <: Rational ? T(string(numerator(val)))/T(string(denominator(val))) : T(string(val))
    elseif T == BigInt
        o = V <: Rational ? T(numerator(val))/T(denominator(val)) : T(val)
    else
        o = convert(T,val)
    end

    return o

end
