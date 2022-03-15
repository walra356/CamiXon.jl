# ======================== sup(i) ==============================================

function _superscript(i::Int)

    c = i < 0 ? [Char(0x207B)] : []

    for d ∈ reverse(digits(abs(i)))
        d == 0 ? push!(c, Char(0x2070)) :
        d == 1 ? push!(c, Char(0x00B9)) :
        d == 2 ? push!(c, Char(0x00B2)) :
        d == 3 ? push!(c, Char(0x00B3)) : push!(c, Char(0x2070+d))
    end

    return join(c)

end

function _subscript(i::Int)

    c = i < 0 ? [Char(0x208B)] : []

    for d ∈ reverse(digits(abs(i)))
        push!(c, Char(0x2080+d))
    end

    return join(c)

end

# ======================== sup(i) ==============================================

@doc raw"""
    sup(i)

Superscript notation for integers and rational numbers
#### Examples:
```
sup(3) * 'P'
 "³P"
```
"""
function sup(i::T) where T<:Real

    sgn = i < 0 ? Char(0x207B) : ""

    num = _superscript(numerator(abs(i)))
    den = _superscript(denominator(abs(i)))

    return T == Rational{Int} ? (sgn * num * '\U141F' * den) : (sgn * num)

end

# ======================== sub(i) ==============================================

@doc raw"""
    sub(i)

Subscript notation for integers and rational numbers
#### Examples:
```
'D' * sub(5//2)
 "D₅⸝₂"
```
"""
function sub(i::T) where T<:Real

    sgn = i < 0 ? Char(0x208B) : ""

    num = _subscript(numerator(abs(i)))
    den = _subscript(denominator(abs(i)))

    return T == Rational{Int} ? (sgn * num * '\U2E1D' * den) : (sgn * num)

end

@doc raw"""
    frac(i)

Fraction notation for rational numbers
#### Examples:
```
frac(-5//2)
 "-⁵/₂"
```
"""
function frac(i::Rational{Int})

    sgn = i < 0 ? "-" : ""

    num = _superscript(numerator(abs(i)))
    den = _subscript(denominator(abs(i)))

    return sgn * num *  '/' * den

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
