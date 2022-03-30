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

# ============================== Value(val, unit) =====================================

"""
    Value(val::Real, unit::String, name::String)

Type with fields:
* ` .val`: numerical value
* `.unit`: unit specifier
* `.name`: type description
* `.symbol`: symbol

Value object
#### Example:
```
f = Value(1,"Hz", "frequency", "f")
 Value(1, "Hz", "frequency")

f.name
 "frequency"
```
"""
struct Value
    val::Number
    unit::String
    name::String
  symbol::String
end

# ========================== strValue(val) =====================================

"""
    strValue(f::Value)

String expression for Value object
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

# ========================= codata =============================================

"""
    Codata

* '∆νCs': Cs hyperfine transition frequency
* '   c': speed of light in vacuum
* '   h': Planck constant
* '   ħ': Planck constant (reduced)
* '   e': elementary charge
* '  kB': Boltzmann constant
* '  NA': Avogadro constant
* ' Kcd': Luminous efficacy
* '  me': electron rest mass
* '  R∞': Rydberg constant
* '  Ry': Rydberg frequency
* '  Eh': Hartree a.u.
* '   α': fine-structure constant
* '  μ0': magnetic permitivity of vacuum
* '  ε0': electric permitivity of vacuum
* '  KJ': Josephson constant
* '  RK': Von Klitzing constant
* '   R': Molar gas constant
* 'matE': unit conversion matrix

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
codata.c.symbol * " = " strValue(codata.c)
 "c = 299792458 m s
```
"""
function createCodata(year::Int)

    year == 2018 || error("Error: codata$(year) not implemented")

    ∆νCs = Value(9192631770, "Hz", sup(133)*"Cs hyperfine transition frequency", "∆νCs")
       c = Value(299792458, "m s"*sup(-1), "speed of light in vacuum", "c")
       h = Value(6.62607015e-34, "J Hz"*sup(-1), "Planck constant", "h")
       ħ = Value(h.val/(2π), "J s", "Planck constant (reduced)", "ħ")
       e = Value(1.602176634e-19, "C", "elementary charge", "e")
      kB = Value(1.380649e-23, "J K"*sup(-1), "Boltzmann constant", "kB")
      NA = Value(6.02214076e23, "mol"*sup(-1), "Avogadro constant", "NA")
     Kcd = Value(683, "lm W"*sup(-1), "Luminous efficacy", "Kcd")
      me = Value(9.1093837015e-31, "Kg", "electron rest mass", "m"*sub2("e"))
      R∞ = Value(10973731.568160, "m"*sup(-1), "Rydberg constant", "R∞")

      Ry = Value(R∞.val*c.val, "Hz", "Rydberg frequency", "Ry")
      Eh = Value(2Ry.val*h.val, "Hartree a.u.", "Hartree atomic unit", "E"*sub2("h"))
       α = Value(sqrt(Eh.val/me.val)/c.val, "", "fine-structure constant", "α")
      μ0 = Value(2α.val*h.val/e.val/e.val/c.val, "N A"*sup(-2), "magnetic permitivity of vacuum", "μ"*sub(0))
      ε0 = Value(1/μ0.val/c.val/c.val, "F m"*sup(-1), "electric permitivity of vacuum", "ε"*sub(0))
      KJ = Value(2e.val/h.val, "Hz V"*sup(-1), "Josephson constant", "KJ")
      RK = Value(h.val/e.val/e.val, "Ω", "Von Klitzing constant", "RK")
       R = Value(NA.val*kB.val, "J mol"*sup(-1)*"K"*sup(-1), "Molar gas constant", "R")

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

    return Codata(∆νCs, c, h, ħ, e, kB, NA, Kcd, me, R∞, Ry, Eh, α, μ0, ε0, KJ, RK, R, M)

end

# =============== convertUnits(val; unitIn="kHz", unitOut="xHz") ===============

"""
    convertUnits(val; unitIn="kHz", unitOut="xHz")

Unit conversion between μHz, ..., EHz, Hartree, Rydberg, Joule, and eV

default input: Hartree

default output: xHz ∈ {μHz, mHz, Hz, kHz, MHz, GHz, THz, PHz, EHz}
#### Example:
```
convertUnits(1; unitIn="Hz", unitOut="Joule")
 6.62607015e-34

convertUnits(1; unitIn="Hartree", unitOut="Hz")
 Value(6.57968392050182e15, "Hz")

f = convertUnits(1) # default input (Hartree) and output (xHz)
strf = strValue(f)
 "6.57968 PHz"
```
"""
function convertUnits(val; unitIn="Hartree", unitOut="xHz", cst=codata)
# ==============================================================================
#  convertUnits(val; unitIn="kHz", unitOut="MHz")
# ==============================================================================
    U = ["μHz","mHz","Hz","kHz","MHz","GHz","THz","PHz","EHz","Hartree","Rydberg","Joule","eV"]

    unitIn ∈ U || error("Error: unitIn ∉ {μHz, ..., EHz, Hartree, Rydberg, Joule, eV}")

    M = cst.matE

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
    calibrationReport(E, Ecal; unitIn="Hartree")

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
function calibrationReport(E, Ecal; unitIn="Hartree")

    T = typeof(E)

    V = BigFloat

    E = myconvert(V, E)
    Ecal = myconvert(V, Ecal)

    ΔE = abs(E-Ecal)
    ΔErel = ΔE/E

    Δf = convertUnits(ΔE)
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
