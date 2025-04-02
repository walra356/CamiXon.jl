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
#                                coadata.jl
# ==============================================================================

# ============================== Value(val, unit) ==============================

@doc raw"""
    Value(val::Real, unit::String)

Object to hold a real numerical value together with a unit specifier.

The fields are:
* ` .val`: numerical value (`::Real`)
* `.unit`: unit specifier (`::String`)
#### Example:
```
julia> f = Value(1,"Hz")
Value(1, "Hz")

julia> f.val
1

julia> f.unit
"Hz"
```
"""
struct Value

    val::Number
    unit::String

end

# ========================== strValue(val) =====================================

@doc raw"""
    strValue(f::Value)

String expression for a [`Value`](@ref) object in `:compact => true`
representation
#### Example:
```
julia> f = Value(1,"Hz")
Value(1, "Hz")

julia> strValue(f)
"1 Hz"
```
"""
function strValue(f::Value)

    strval = repr(f.val, context=:compact => true)
    strunit = " " * f.unit

    return strval * strunit

end

# ============================== NamedValue(val, unit) =========================

@doc raw"""
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
julia> f = Value(1,"Hz")
Value(1, "Hz")

julia> f = castNamedValue(f, name="frequency", comment="comment")
NamedValue(Value(1, "Hz"), "frequency", "comment")

julia> f.name
"frequency"
```
"""
struct NamedValue

    val::Value
    name::String
    comment::String

end

# ==================== NamedValue(val, unit) ===================================

@doc raw"""
    castNamedValue(val::Value; name=" ", comment=" ")

Method to create a [`NamedValue`](@ref) object
#### Example
```
julia> v = Value(1.602176634e-19, "C");

julia> nv = castNamedValue(v; name="e");

julia> nv.name * " = " * strValue(nv.val)
"e = 1.60218e-19 C"
```
"""
function castNamedValue(val::Value; name=" ", comment=" ")

    return NamedValue(val, name, comment)

end

# ========================= codata =============================================

@doc raw"""
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
* `   .u`: unified atomic mass unit (`::Value`)
* `.matE`: unit conversion matrix (Matrix{Float64})

#### Example:
```
julia> codata = castCodata(2022);

julia> codata.μ0
Value(1.2566370612696005e-6, "N A⁻²")

julia> codata.μ0.val
1.2566370612696005e-6
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
      mp::Value
      R∞::Value
      Ry::Value
      Eh::Value
       α::Value
      a0::Value
      μB::Value
      μN::Value
      μ0::Value
      ε0::Value
      KJ::Value
      RK::Value
       R::Value
       u::Value
    matE::Matrix{Float64}

end

# ========================= castCodata(year) =================================

@doc raw"""
    castCodata(year::Int)

Method to create the [`Codata`](@ref) object
#### Example:
```
julia> codata = castCodata(2022);

julia> strValue.([codata.∆νCs,codata.c,codata.h])
3-element Vector{String}:
 "9192631770 Hz"
 "299792458 m s⁻¹"
 "6.62607e-34 J Hz⁻¹"
```
"""
function castCodata(year::Int)

    if year == 2022
        ∆νCs = Value(9192631770, "Hz")
        c = Value(299792458, "m s" * sup(-1))
        h = Value(6.62607015e-34, "J Hz" * sup(-1))
        ħ = Value(h.val / (2π), "J s")
        e = Value(1.602176634e-19, "C")
        kB = Value(1.380649e-23, "J K" * sup(-1))
        NA = Value(6.02214076e23, "mol" * sup(-1))
        Kcd = Value(683, "lm W" * sup(-1))
        me = Value(9.1093837139e-31, "kg")
        mp = Value(1.67262192595e-27, "kg")
        R∞ = Value(10973731.568157, "m" * sup(-1))
        u = Value(1.66053906892e-27, "kg")
    elseif year == 2018
        ∆νCs = Value(9192631770, "Hz")
        c = Value(299792458, "m s" * sup(-1))
        h = Value(6.62607015e-34, "J Hz" * sup(-1))
        ħ = Value(h.val / (2π), "J s")
        e = Value(1.602176634e-19, "C")
        kB = Value(1.380649e-23, "J K" * sup(-1))
        NA = Value(6.02214076e23, "mol" * sup(-1))
        Kcd = Value(683, "lm W" * sup(-1))
        me = Value(9.1093837015e-31, "kg")
        mp = Value(1.67262192595e-27, "kg")
        R∞ = Value(10973731.568160, "m" * sup(-1))
        u = Value(1.66053906660e-27, "kg")
    else
       #error("Error: codata$(year) not implemented - use codata 2018 or 2022 ")
        throw(DomainError(year))
    end

    Ry = Value(R∞.val * c.val, "Hz")
    Eh = Value(2Ry.val * h.val, "J")
    α = Value(sqrt(Eh.val / me.val) / c.val, "")
    μ0 = Value(2α.val * h.val / e.val / e.val / c.val, "N A" * sup(-2))
    ε0 = Value(1 / μ0.val / c.val / c.val, "F m" * sup(-1))
    KJ = Value(2e.val / h.val, "Hz V" * sup(-1))
    RK = Value(h.val / e.val / e.val, "Ω")
    R = Value(NA.val * kB.val, "J mol" * sup(-1) * "K" * sup(-1))
    a0 = Value(ħ.val/ α.val / me.val/ c.val, "m")
    μB = Value(e.val * ħ.val / 2.0 / me.val, "J T" * sup(-1))
    μN = Value(e.val * ħ.val / 2.0 / mp.val, "J T" * sup(-1))

    H = 1e-18 * c.val * 2R∞.val
    J = 1e+18 * h.val
    V = 1e+18 * h.val / e.val
    K = h.val/kB.val
    Ehv = Eh.val
    kBv = kB.val
    ev = e.val
    cv = c.val
    hc = h.val * c.val
    

    M = [
#   μHz mHz Hz  kHz MHz  GHz  THz  PHz  EHz  Eh     Ry     J      eV      cm-1      m-1     K          mK        μK
    1   1e3 1e6 1e9 1e12 1e15 1e18 1e21 1e24 1e24*H 2e23*H 1e24/J 1e24/V  1e8cv     1e10cv  1e-6/K     1e-9/K    1e-12/K;       # 1 - μHz
    1   1   1e3 1e6 1e9  1e12 1e15 1e18 1e21 1e21*H 2e20*H 1e21/J 1e21/V  1e5cv     1e7cv   1e-3/K     1e-6/K    1e-9/K;        # 2 - mHz
    1   1   1   1e3 1e6  1e9  1e12 1e15 1e18 1e18*H 2e17*H 1e18/J 1e18/V  1e2cv     1e5cv   1/K        1e-3/K    1e-6/K;        # 3 - Hz
    1   1   1   1   1e3  1e6  1e9  1e12 1e15 1e15*H 2e14*H 1e15/J 1e15/V  1e-1cv    1e1cv   1e3/K      1/K       1e-3/K;        # 4 - kHz
    1   1   1   1   1    1e3  1e6  1e9  1e12 1e12*H 2e11*H 1e12/J 1e12/V  1e-4cv    1e-2cv  1e6/K      1e3/K     1/K;           # 5 - MHz
    1   1   1   1   1    1    1e3  1e6  1e9  1e9*H  2e8*H  1e9/J  1e9/V   1e-7cv    1e-5cv  1e9/K      1e6/K     1e3/K;         # 6 - GHz
    1   1   1   1   1    1    1    1e3  1e6  1e6*H  2e5*H  1e6/J  1e6/V   1e-10cv   1e-8cv  1e12/K     1e9/K     1e6/K;         # 7 - THz
    1   1   1   1   1    1    1    1    1e3  1e3*H  2e2*H  1e3/J  1e3/V   1.e-13cv  1e-11cv 1e15/K     1e12/K    1e9/K;         # 8 - PHz
    1   1   1   1   1    1    1    1    1    H      0.5*H  1/J    1/V     1.e-16cv  1e-14cv 1e18/K     1e15/K    1e12/K;        # 9 - EHz
    1   1   1   1   1    1    1    1    1    1      0.5    1/Ehv  ev/Ehv  100hc/Ehv hc/Ehv  kBv/Ehv    1e-3kBv/Ehv 1e-6kBv/Ehv; # 10 - Hartree
    1   1   1   1   1    1    1    1    1    1      1      2/Ehv  2ev/Ehv 200hc/Ehv 2hc/Ehv 2kBv/Ehv   2e-3kBv/Ehv 2e-6kBv/Ehv; # 11 - Rydberg
    1   1   1   1   1    1    1    1    1    1      1      1      e.val   100hc     hc      kBv        1e-3kBv     1e-6kBv;     # 12 - Joule
    1   1   1   1   1    1    1    1    1    1      1      1      1       100hc/ev  hc/ev   kBv/ev     1e-3kBv/ev  1e-6kBv/ev;  # 13 - eV
    1   1   1   1   1    1    1    1    1    1      1      1      1       1         0.01    0.01kBv/hc 1e-5kBv/hc  1e-8kBv/hc ; # 14 cm-1
    1   1   1   1   1    1    1    1    1    1      1      1      1       1         1       kBv/hc     1e-3kBv/hc  1e-6kBv/hc;  # 15 m-1
    1   1   1   1   1    1    1    1    1    1      1      1      1       1         1          1       1e-3        1e-6;        # 16 K
    1   1   1   1   1    1    1    1    1    1      1      1      1       1         1          1       1           1e-3;        # 17 mK
    1   1   1   1   1    1    1    1    1    1      1      1      1       1         1          1       1           1            # 18 μK
 
    ]


    N = Int(sqrt(length(M)))

    for i=1:N
        for j=i:N
            M[j,i] = 1/M[i,j]
        end
    end

    o = Codata(∆νCs,c,h,ħ,e,kB,NA,Kcd,me,mp,R∞,Ry,Eh,α,a0,μB,μN,μ0,ε0,KJ,RK,R,u,M)

    return o

end

# ========================= listCodata(codata::Codata) =========================

@doc raw"""
    listCodata(codata::Codata)

Method to list the fields of [`Codata`](@ref) by their symbolic name
#### Example:
```
julia> julia> codata = castCodata(2018);

julia> listCodata(codata)
∆νCs = 9192631770 Hz      - ¹³³Cs hyperfine transition frequency
   c = 299792458 m s⁻¹    - speed of light in vacuum
   h = 6.62607e-34 J Hz⁻¹ - Planck constant
   ħ = 1.05457e-34 J s    - Planck constant (reduced)
   e = 1.60218e-19 C      - elementary charge
  kB = 1.38065e-23 J K⁻¹  - Boltzmann constant
  NA = 6.02214e23 mol⁻¹   - Avogadro constant
 Kcd = 683 lm W⁻¹         - Luminous efficacy
  mₑ = 9.10938e-31 kg     - electron mass
  mₚ = 1.67262e-27 kg     - proton mass
  R∞ = 1.09737e7 m⁻¹      - Rydberg constant
  Ry = 3.28984e15 Hz      - Rydberg frequency
  Eₕ = 4.35974e-18 J      - Hartree atomic unit
   α = 0.00729735         - fine-structure constant
  a0 = 5.29177e-11 m      - Bohr radius
  μB = 9.27401e-24 J T⁻¹  - Bohr magneton
  μN = 5.05078e-27 J T⁻¹  - nuclear magneton
  μ₀ = 1.25664e-6 N A⁻²   - magnetic permitivity of vacuum
  ε₀ = 8.85419e-12 F m⁻¹  - electric permitivity of vacuum
  KJ = 4.83598e14 Hz V⁻¹  - Josephson constant
  RK = 25812.8 Ω          - Von Klitzing constant
   R = 8.31446 J mol⁻¹K⁻¹ - Molar gas constant
   u = 1.66054e-27 kg     - unified atomic mass unit
```
"""
function listCodata(codata::Codata; msg=true)

    ∆νCs = NamedValue(codata.∆νCs, "∆νCs",
                       sup(133)*"Cs hyperfine transition frequency")
       c = NamedValue(codata.c, "c", "speed of light in vacuum")
       h = NamedValue(codata.h, "h", "Planck constant")
       ħ = NamedValue(codata.ħ, "ħ", "Planck constant (reduced)")
       e = NamedValue(codata.e, "e", "elementary charge")
      kB = NamedValue(codata.kB, "kB", "Boltzmann constant")
      NA = NamedValue(codata.NA, "NA", "Avogadro constant")
     Kcd = NamedValue(codata.Kcd, "Kcd", "Luminous efficacy")
      me = NamedValue(codata.me, "m"*sub("e"), "electron mass")
      mp = NamedValue(codata.mp, "m"*sub("p"), "proton mass")
      R∞ = NamedValue(codata.R∞, "R∞", "Rydberg constant")
      Ry = NamedValue(codata.Ry, "Ry", "Rydberg frequency")
      Eh = NamedValue(codata.Eh, "E"*sub("h"), "Hartree atomic unit")
       α = NamedValue(codata.α, "α", "fine-structure constant")
      a0 = NamedValue(codata.a0, "a0", "Bohr radius")
      μB = NamedValue(codata.μB, "μB", "Bohr magneton")
      μN = NamedValue(codata.μN, "μN", "nuclear magneton")
      μ0 = NamedValue(codata.μ0, "μ"*sub(0), "magnetic permitivity of vacuum")
      ε0 = NamedValue(codata.ε0, "ε"*sub(0), "electric permitivity of vacuum")
      KJ = NamedValue(codata.KJ, "KJ", "Josephson constant")
      RK = NamedValue(codata.RK, "RK", "Von Klitzing constant")
       R = NamedValue(codata.R, "R", "Molar gas constant")
       u = NamedValue(codata.u, "u", "unified atomic mass unit")

    U =  [∆νCs, c, h, ħ, e, kB, NA, Kcd, me, mp, R∞, Ry, Eh, α, a0, μB, μN, μ0, ε0, KJ, RK, R, u]

    str = ""
    for i ∈ eachindex(U)
        str *= lpad(U[i].name, 4) * " = " * rpad(strValue(U[i].val), 18) * " - " * U[i].comment * "\n"
    end

    return msg ? println(str) : str

end

# ===== convertUnit(val, codata; unitIn="kHz", unitOut="xHz") ===========

@doc raw"""
    convertUnit(val, codata; unitIn="Hartree", unitOut="xHz")

Unit conversion between μHz,⋯ EHz, Hartree, Rydberg, J eV, cm-1, m-1, K, mK μK

default input: Hartree

default output: xHz ∈ {μHz, mHz, Hz, kHz, MHz, GHz, THz, PHz, EHz}
#### Example:
```
julia> codata = castCodata(2018);
julia> convertUnit(1, codata; unitIn="Hz", unitOut="J")
Value(6.62607015e-34, "J")

julia> convertUnit(1, codata; unitIn="Hartree", unitOut="Hz")
Value(6.579683920501762e15, "Hz")

julia> f = convertUnit(1, codata) # default input (Hartree) and output (xHz);
julia> strf = strValue(f)
"6.57968 PHz"
```
"""
function convertUnit(val, codata; unitIn="Hartree", unitOut="xHz")
    # ==============================================================================
    #  convertUnit(val, codata; unitIn="kHz", unitOut="MHz")
    # ==============================================================================
        U = ["μHz","mHz","Hz","kHz","MHz","GHz","THz","PHz","EHz"]
            Base.push!(U,"Hartree")
            Base.push!(U,"Rydberg")
            Base.push!(U,"J")
            Base.push!(U,"eV")
            Base.push!(U,"cm-1")
            Base.push!(U,"m-1")
            Base.push!(U,"K")
            Base.push!(U,"mK")
            Base.push!(U,"μK")
    
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

@doc raw"""
    calibrationReport(E, Ecal, codata::Codata; unitIn="Hartree", msg=true)

Comparison of energy E with calibration value Ecal

default input: Hartree
#### Example:
```
julia> codata = castCodata(2022)
julia> calibrationReport(1.1, 1.0, codata; unitIn="Hartree")

calibration report (Float64):
Ecal = 1 Hartree
E = 1.1000000000000001 Hartree
absolute accuracy: ΔE = 0.1 Hartree (657.968 THz)
relative accuracy: ΔE/E = 0.0909091
```
"""
function calibrationReport(E, Ecal, codata::Codata; unitIn="Hartree", msg=true)

    T = typeof(E)

    B = BigFloat

    E = convert(B, E)
    Ecal = convert(B, Ecal)

    ΔE = abs(E-Ecal)
    ΔErel = ΔE/E

    Δf = convertUnit(ΔE, codata)
    strΔf = " (" * strValue(Δf) * ")"
    strΔE = repr(ΔE, context=:compact => true)
    strΔErel = repr(ΔErel, context=:compact => true)

    str = "\ncalibration report ($T):\n"
    str *= Printf.@sprintf "Ecal = %.17g %s \n" Ecal unitIn
    T == B ? (str *= Printf.@sprintf "E = %.25g %s \n" E unitIn) : (str *= Printf.@sprintf "E = %.17g %s \n" E unitIn)
    str *= ΔE ≠ 0 ? "absolute accuracy: ΔE = " * strΔE * " " * unitIn * strΔf * "\n" :
                    "absolute accuracy: ΔE = 0 (exact under $T precision)\n"
    str *= ΔE ≠ 0 ? "relative accuracy: ΔE/E = " * strΔErel * "\n"                   :
                    "relative accuracy: ΔE/E = 0 (exact under $T precision)\n"

    return msg ? println(str) : str

end
