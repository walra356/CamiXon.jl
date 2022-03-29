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

# ============================== Value(val, unit) =====================================

"""
    Value(val::Real, unit::String)

Type with fields:
* ` .val`: numerical value
* `.unit`: unit specifier

Value object
#### Example:
```
f = Value(1,"Hz")
str = "$(f.val) $(f.unit)"
 "1 Hz"
```
"""
struct Value
    val::Number
    unit::String
end

# ========================== strValue(val) ===================================

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
function convertUnits(val; unitIn="Hartree", unitOut="xHz")
# ==============================================================================
#  convertUnits(val; unitIn="kHz", unitOut="MHz")
# ==============================================================================
    U = ["μHz","mHz","Hz","kHz","MHz","GHz","THz","PHz","EHz","Hartree","Rydberg","Joule","eV"]

    unitIn ∈ U || error("Error: unitIn ∉ {μHz, ..., EHz, Hartree, Rydberg, Joule, eV}")

    N = length(U)

    H = 0.15198298460570
    J = 6.62607015e-19
    V = 4.135667696

    M =[1 1e3 1e6 1e9 1e12 1e15 1e18 1e21 1e24 1e+21/H 0.5e+21/H 1e21/J 1e21/V;
        1 1 1e3 1e6 1e9 1e12 1e15 1e18 1e21 1e+18/H 0.5e+18/H 1e18/J 1e18/V;
        1 1 1 1e3 1e6 1e9 1e12 1e15 1e18 1e+15/H 0.5e+15/H 1e15/J 1e15/V;
        1 1 1 1 1e3 1e6 1e9 1e12 1e15 1e+12/H 0.5e+12/H 1e12/J 1e12/V;
        1 1 1 1 1 1e3 1e6 1e9 1e12 1e+9/H 0.5e+9/H 1e9/J 1e9/V;
        1 1 1 1 1 1 1e3 1e6 1e9 1e+6/H 0.5e+6/H 1e6/J 1e6/V;
        1 1 1 1 1 1 1 1e3 1e6 1e3/H 0.5e+3/H 1e3/J 1e3/V;
        1 1 1 1 1 1 1 1 1e3 1/H 0.5/H 1/J 1/V;
        1 1 1 1 1 1 1 1 1 1e-3/H 0.5e-3/H 1e-3/J 1e-3/V;
        1 1 1 1 1 1 1 1 1 1 0.5 2.2937122783963e17 1/27.211386245988;
        1 1 1 1 1 1 1 1 1 1 1 2*2.2937122783963e17 2/27.211386245988;
        1 1 1 1 1 1 1 1 1 1 1 1 1/6.241509074e18;
        1 1 1 1 1 1 1 1 1 1 1 1 1]

    for i=1:N
        for j=i:N
            M[j,i] = 1/M[i,j]
        end
    end

    v = zeros(Float64,N)
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
