"""
    find_all(A [,a...]; count=false)
A: string/array of elements of the same type

default   : Array containing the index (indices) of selected elements of A (default: all elements)

count=true: The number of indices found for selected elements of A (default: all elements)
#### Examples:
```
A = [:📑,:📌,:📢,:📌,:📞]
B = [1,2,3,2,5]
str = "aβcβd";
find_all(A) == find_all(B) == find_all(str)
true

find_all(A,:📌)
1-element Array{Array{Int64,1},1}:
 [2, 4]

find_all(str)
4-element Array{Array{Int64,1},1}:
 [1]
 [2, 4]
 [3]
 [5]

find_all(A; count=true)
4-element Array{Int64,1}:
 1
 2
 1
 1

str = "📑📌📢📌📞"
find_all(str,'📌')
1-element Array{Array{Int64,1},1}:
 [2, 4]
```
"""
function find_all(A::Union{String,AbstractArray{T,1}}, a::T...; count=false)  where T

    typeof(A) == String ? A = Base.collect(A) : nothing

    a == () ? a = Base.unique(A) : nothing

    o = [Base.findall(x -> x == a[i], A) for i in eachindex(a)]

    return count ? length.(o) : o

end


"""
    find_first(A [,a...]; dict=false)
The first index of selected Array element

A: string/array of elements of the same type

default  : Array containing the first index (indices) of selected elements of A (default: all elements)

dict=true: Dict for the first index (indices) of selected elements of A (default: all elements)
#### Examples:
```
A = [:📑,:📌,:📢,:📌,:📞]
B = [1,2,3,2,5]
str = "aβcβd";

find_first(A) == find_first(B) == find_first(str)
true

find_first(A,:📌)
1-element Array{Array{Int64,1},1}:
 2

find_last(A,:📌; dict=true)
1-element Array{Pair{Symbol,Int64},1}:
 :📌 => 2

find_last(A; dict=true)
4-element Array{Pair{Symbol,Int64},1}:
 :📑 => 1
 :📌 => 2
 :📢 => 3
 :📞 => 5

find_first(str)
4-element Array{Int64,1}:
 1
 2
 3
 5
```
"""
function find_first(A::Union{String,AbstractArray{T,1}}, a::T...; dict=false)  where T

    typeof(A) == String ? A = Base.collect(A) : nothing

    a == () ? a = Base.unique(A) : nothing

    o = [Base.findfirst(x -> x == a[i], A) for i in eachindex(a)]

    return dict ? [a[i] => o[i] for i in eachindex(a)] : o

end

"""
    find_last(A [,a...]; dict=false)

The last index of selected Array element

A: string/array of elements of the same type

default  : Array containing the lasst index (indices) of selected elements of A (default: all elements)

dict=true: Dict for the lasst index (indices) of selected elements of A (default: all elements)
#### Examples:
```
A = [:📑,:📌,:📢,:📌,:📞]
B = [1,2,3,2,5]
str = "aβcβd";
find_last(A) == find_first(B) == find_first(str)
true

find_last(A,:📌)
1-element Array{Array{Int64,1},1}:
 4

find_last(A,:📌; dict=true)
1-element Array{Pair{Symbol,Int64},1}:
 :📌 => 4

find_last(A; dict=true)
4-element Array{Pair{Symbol,Int64},1}:
 :📑 => 1
 :📌 => 4
 :📢 => 3
 :📞 => 5

find_last(str)
4-element Array{Int64,1}:
 1
 4
 3
 5
```
"""
function find_last(A::Union{String,AbstractArray{T,1}}, a::T...; dict=false)  where T

    typeof(A) == String ? A = Base.collect(A) : nothing

    a == () ? a = Base.unique(A) : nothing

    o = [Base.findlast(x -> x == a[i], A) for i ∈ eachindex(a)]

    return dict ? [a[i] => o[i] for i ∈ eachindex(a)] : o

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

# =============== convertUnits(val; unitIn="kHz", unitOut="MHz") ===============

"""
    convertUnits(val; unitIn="kHz", unitOut="MHz")

Unit conversion between μHz, mHz, Hz, kHz, MHz, GHz, THz, PHz, EHz, Hartree, Rydberg, Joule, and eV
see also frequencyUnits(val; unitIn="Hartree")
#### Example:
```
convertUnits(1; unitIn="Hz", unitOut="Joule")
 6.62607015e-34
```
"""
function convertUnits(val; unitIn="Hartree", unitOut="Hz")
# ==============================================================================
#  convertUnits(val; unitIn="kHz", unitOut="MHz")
# ==============================================================================
    U = ["μHz","mHz","Hz","kHz","MHz","GHz","THz","PHz","EHz","Hartree","Rydberg","Joule","eV"]

    unitIn  ∈ U || error("Error: conversion not implemented")
    unitOut ∈ U || error("Error: conversion not implemented")

    N = length(U)
    H = 0.15198298460570
    J = 6.62607015e-19
    V = 4.135667696

    M =[1 1e3 1e6 1e9 1e12 1e15 1e18 1e21 1e24 1e+21/H 2e+21/H 1e21/J 1e21/V;
        1 1 1e3 1e6 1e9 1e12 1e15 1e18 1e21 1e+18/H 2e+18/H 1e18/J 1e18/V;
        1 1 1 1e3 1e6 1e9 1e12 1e15 1e18 1e+15/H 2e+15/H 1e15/J 1e15/V;
        1 1 1 1 1e3 1e6 1e9 1e12 1e15 1e+12/H 2e+12/H 1e12/J 1e12/V;
        1 1 1 1 1 1e3 1e6 1e9 1e12 1e+9/H 2e+9/H 1e9/J 1e9/V;
        1 1 1 1 1 1 1e3 1e6 1e9 1e+6/H 2e+6/H 1e6/J 1e6/V;
        1 1 1 1 1 1 1 1e3 1e6 1e3/H 2e+3/H 1e3/J 1e3/V;
        1 1 1 1 1 1 1 1 1e3 1/H 2/H 1/J 1/V;
        1 1 1 1 1 1 1 1 1 1e-3/H 2e-3/H 1e-3/J 1e-3/V;
        1 1 1 1 1 1 1 1 1 1 0.5 2.2937122783963e17 1/27.211386245988;
        1 1 1 1 1 1 1 1 1 1 1 2*2.2937122783963e17 2/27.211386245988;
        1 1 1 1 1 1 1 1 1 1 1 1 1/6.241509074e18;
        1 1 1 1 1 1 1 1 1 1 1 1 1]

    for i=1:N
        for j=i:N
           M[j,i] = 1/M[i,j]
        end
    end

    id = findfirst(x -> x==unitIn, U)

    v = zeros(Float64,N)
    v[id] = val

    w = M * v

    id = findfirst(x -> x==unitOut, U)

    return w[id]

end

# ============================== Frequency =====================================

"""
    Frequency(val::Real, unit::String)

Type with fields:
* ` .val`: frequency value
* `.unit`: frequency unit

Frequency object
#### Example:
```
f = Frequency(1, "Hz")
f.val
 1

f.unit
 "Hz"
```
"""
struct Frequency
    val::Real
    unit::String
end

# ========================== strFrequency(f) ===================================

"""
    strFrequency(f::Frequency)

String for frequency object
#### Example:
```
f = Frequency(1, "Hz")
strFrequency(f)
 "1 Hz"
```
"""
function strFrequency(f::Frequency)

    strval = repr(f.val, context=:compact => true)
    strunit = " " * f.unit

    return strval * strunit

end

# ========================== frequencyUnits(f) =================================

"""
    frequencyUnits(val; unitIn="Hartree")

Energy in frequency units - see also: convertUnits(val; unitIn="Hartree", unitOut="Hz")
#### Example:
```
f = frequencyUnits(1; unitIn="Hartree")
 Frequency(6.57968392050182, "PHz")

strFrequency(f)
 "6.57968 PHz"
```
"""
function frequencyUnits(val; unitIn="Hartree")

    U = ["Hartree"]

    unitIn ∈ U || error("Error: conversion not implemented")

    mul = 0.6579683920502001

    unitOut = mul * 1e3 < val ? "EHz" :
        mul * 1e-0  ≤ val < mul * 1e3   ? "PHz" :
        mul * 1e-3  ≤ val < mul ? "THz" :
        mul * 1e-7  ≤ val < mul * 1e-3  ? "GHz" :
        mul * 1e-10 ≤ val < mul * 1e-7  ? "MHz" :
        mul * 1e-13 ≤ val < mul * 1e-10 ? "kHz" :
        mul * 1e-16 ≤ val < mul * 1e-13 ? "Hz"  :
        mul * 1e-19 ≤ val < mul * 1e-16 ? "mHz" : "μHz"

    f = convertUnits(val; unitIn, unitOut)

    return Frequency(f, unitOut)

end

# ============ calibrationReport(E, Ecal; unitIn="Hartree") ====================

"""
    calibrationReport(E, Ecal; unitIn="Hartree")

Comparison of energy E with calibration value Ecal
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

    Δf = frequencyUnits(ΔE; unitIn="Hartree")
    strΔf = strFrequency(Δf)
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
