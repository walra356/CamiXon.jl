# ==================================== polynomial(c,x;deriv) ===================

@doc raw"""
    polynomial(coords::Vector{T}, x::T[; deriv=0]) where T<:Number

Method to evaluate the function ``f(x)=\text{polynomial}(c,x)``, where
``c=[c_0,⋯\ c_d]`` is the vector representation of a polynomial of degree ``d``.
```math
    \text{polynomial}(c,x)=c_0 + c_1 x + ⋯ + c_d x^d.
```
### Examples:
```
coords = ones(Int,6)                     # for polynomial of degree 5 with unit coefficients
f0(x) = polynomial(coords,x)             # default
fd(x) = polynomial(coords,x; deriv=1)    # first derivative
fp(x) = polynomial(coords,x; deriv=-1)   # primitive (with zero integration constant)
f0(1)
 6

fd(1)
 15

fp(1)
 49//20
```
"""
function polynomial(coords::Vector{T}, x::V; deriv=0) where {T <: Real, V <: Real}

    coords = deriv == 0 ? coords : deriv ≥ Base.length(coords) ? 0 : deriv == -1 ? polynom_primitive(coords) : polynom_derivatives(coords; deriv)

    k = Base.length(coords)
    X = Base.ones(T,k)

    for i=2:k
        X[i] = X[i-1] * x
    end

    o = LinearAlgebra.dot(coords, X)

    return o

end

# ==================================== polynom_derivative(coords) ==============

@doc raw"""
    polynom_derivative(coords)

Vector representation of the first derivative of the polynomial `coords`,
```math
    p'(c,x)=c_1 + 2 c_2 x + ⋯ + d c_d x^{d-1},
```
Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``.
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,⋯\ c_d]``
consisting of the polynomial coefficients.
### Examples:
```
coords=[1,1,1,1,1]                 # vector representation of polynomial of degree d=4
polynom_derivative(coords)         # (first) derivative of polynomial `coords`
4-element Vector{Int64}:
 1
 2
 3
 4
```
"""
function polynom_derivative(coords::Vector{T}) where T<:Real

    k = Base.length(coords)
    k > 1 || return [0]

    return coords[2:end] .* Base.OneTo(k-1)

end

# =============================== polynom_derivative(coords[,deriv=0]) =========

@doc raw"""
    polynom_derivatives(coords::Vector{T}; deriv=0) where T<:Real

Vector representation of derivatives of the polynomial `coords`.

Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``.
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,⋯\ c_d]``
consisting of the polynomial coefficients.

`deriv`: derivative of choice; `default`: `coords` remains unchanged.

### Examples:
```
coords=[1,1,1,1,1]               # vector representation of a polynomial of degree d=4
polynom_derivatives(coords)      # default no (zero) derivative of polynomial `coords`
5-element Vector{Vector{Int64}}:
 1
 1
 1
 1
 1

polynom_derivatives(coords; deriv=2)        # second derivative of polynomial `coords`
3-element Vector{Int64}:
  2
  6
 12
```
"""
function polynom_derivatives(coords::Vector{T}; deriv=0) where T<:Real

    deriv < 0 && error("jwError: negative derivative not defined")

    for k=1:deriv
        coords = CamiXon.polynom_derivative(coords)
    end

    return coords

end

# =============================== polynom_derivative_all(coords) =========

@doc raw"""
    polynom_derivatives_all(coords::Vector{<:Number})

Vector representation of all nontrivial derivatives of the polynomial `coords`.

Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``.
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,⋯\ c_d]``
consisting of the polynomial coefficients.

### Examples:
```
coords=[1,1,1,1,1]               # vector representation of a polynomial of degree d=4
polynom_derivatives_all(coords)      # `all' (nontrivial) derivatives of polynomial `coords`
5-element Vector{Vector{Int64}}:
 [1, 2, 3, 4]
 [2, 6, 12]
 [6, 24]
 [24]
```
"""
function polynom_derivatives_all(coords::Vector{T}) where T<:Real

    k = Base.length(coords)

    coords = CamiXon.polynom_derivative(coords)

    o = [coords]

    for i=2:k-1
        coords = polynom_derivative(coords)
        Base.push!(o,coords)
    end

    return o

end
# ==================================== polynom_power(coords, p) ================

@doc raw"""
    polynom_power(coords, p)

Vector representation of the polynomial `coords` raised to the power `p` which
results in a polynomial in a vector space of dimension ``p d + 1``.

Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``.
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,⋯\ c_d]``
consisting of the polynomial coefficients.
### Examples:
```
coords=[1,1,1]             # vector representation of polynomial of degree ``d=2``
polynom_power(coords,2)
5-element Vector{Int64}:
 1
 2
 3
 2
 1
```
"""
function polynom_power(coords::Vector{T}, power::Int) where T<:Real

    power >= 0 || error("jwError: negative powers not allowed")
    power == 2 && return polynom_product(coords, coords)
    power == 1 && return coords
    power == 0 && return [1]

    o = CamiXon.polynom_product(coords, coords)

    for i=1:power-2
        o = CamiXon.polynom_product(o, coords)
    end

    return o

end

# ==================================== polynom_powers(coords, pmax) ============

@doc raw"""
    polynom_powers(coords::Vector{T}, pmax::Int) where T<:Real

The polynomial `coords` raised to the powers 1,...,pmax  which
results in a collection of polynomials in vector spaces of dimension ``d+1`` tot ``p d + 1``.

Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``.
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,⋯\ c_d]``
consisting of the polynomial coefficients.
### Examples:
```
coords=[1,1,1]                   # vector representation of polynomial of degree d=2
polynom_powers(coords,3)
3-element Vector{Vector{Int64}}:
 [1, 1, 1]
 [1, 2, 3, 2, 1]
 [1, 3, 6, 7, 6, 3, 1]
```
"""
function polynom_powers(coords::Vector{T}, pmax::Int) where T<:Real

    pmax > 0 || error("jwError: minimum power included is unity")

    o = [coords]

    for i=1:pmax-1
        Base.push!(o,CamiXon.polynom_product(o[end], coords))
    end

    return o

end

# ==================================== polynom_primitive(coords) ====================

@doc raw"""
    polynom_primitive(coords::Vector{T}) where T<:Real

Vector representation of the primitive of the polynomial `coords` which is a
polynomial in a vector space of dimension ``p d + 1``.
```math
    P(c,x)=c_{int} +c_0 x + \frac{1}{2} c_1 x^2 + \frac{1}{3} c_2 x^3 + ⋯ + \frac{1}{d+1} c_d x^{d+1},
```
The constant of integration is set to zero, ``c_{int} = 0``.

Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``.
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,⋯\ c_d]``
consisting of the polynomial coefficients.
### Examples:
```
coords=[1,1,1,1,1]         # vector representation of polynomial of degree ``d=4``
polynom_primitive(coords)
6-element Vector{Rational{Int64}}:
 0//1
 1//1
 1//2
 1//3
 1//4
 1//5
```
"""
function polynom_primitive(coords::Vector{T}) where T<:Real

    d = [1//p for p ∈ Base.eachindex(coords)]

    coords = coords .* d

    return Base.pushfirst!(coords,0)      # constant of integration equal to zero

end

# ==================================== polynom_product(a, b) ============================================================

@doc raw"""
    polynom_product(a::Vector{T}, b::Vector{V}) where {T<:Real, V<:Real}

Vector representation of the product of two polynomials, ``a`` and ``b`` which
is a polynomial in a vector space of dimension ``d=m+n``,
```math
    p(c,x)=a_0b_0 + (a_0b_1 + b_0a_1)x + ⋯ + a_n b_m x^{n+m}.
```
Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,⋯\ c_d]``
consisting of the polynomial coefficients.
####
```
[polynom_product1([1.0,1],[1,-1,2])]
 [1.0, 0.0, 1.0, 2.0]

[polynom_product1([1//1,1],[1,-1,2])]
 [1//1, 0//1, 1//1, 2//1]

[polynom_product([1,1],[1,- 1,2])]
 [1, 0, 1, 2]

[polynom_product([1,- 1,2],[1,1])]
 [1, 0, 1, 2]
```
"""
function polynom_product(a::Vector{T}, b::Vector{V}) where {T<:Real, V<:Real}

    n = Base.length(a)
    m = Base.length(b)

    a,b = Base.promote(a,b)

    if m ≥ n
        o = [Base.sum(a[1+j-i]*b[1+i] for i=0:j) for j=0:n-1]
        if m≠n Base.append!(o,[Base.sum(a[n-i]*b[1+i+j] for i=0:n-1) for j=1:m-n]) end
        Base.append!(o,[Base.sum(a[n-i]*b[1+i+j+m-n] for i=0:n-1-j) for j=1:n-1])
    else
        o = [Base.sum(b[1+j-i]*a[1+i] for i=0:j) for j=0:m-1]
        if m≠n Base.append!(o,[Base.sum(b[m-i]*a[1+i+j] for i=0:m-1) for j=1:n-m]) end
        Base.append!(o,[Base.sum(b[m-i]*a[1+i+j+n-m] for i=0:m-1-j) for j=1:m-1])
    end

    return o

end

# ==================================== polynom_product_expansion(a, b, p) ============================================================

@doc raw"""
    polynom_product_expansion(a::Vector{T}, b::Vector{T}, p::Int) where T<:Real

Vector representation of the product of two polynomials, ``a`` (of degree ``n``) and ``b`` (of degree ``m``), with ``m≤n``
truncated at the order ``p`` is a polynomial in a vector space of dimension ``d=p+1``. If ``ab`` is the `polynom_product`,
the `polynom_product_expansion` is ``ab[1:p+1]``
####
```
a = [1,-1,1]
b = [1,1,-1,1,1,1]
o = polynom_product(a, b); println(o)
 [1, 0, -1, 3, -1, 1, 0, 1]

o = expand_product(a, b, 4); println(o)
 [1, 0, -1, 3, -1]

```
"""
function polynom_product_expansion(a::Vector{T}, b::Vector{T}, p::Int) where T<:Real

    n = Base.length(a)
    m = Base.length(b)

    if m ≥ n
        o = [Base.sum(a[1+j-i]*b[1+i] for i=0:j) for j=0:min(n-1,p)]
        p+1 == length(o) && return o
        if m≠n Base.append!(o,[Base.sum(a[n-i]*b[1+i+j] for i=0:n-1) for j=1:min(m-n,p-n+1)]) end
        p+1 == length(o) && return o
        Base.append!(o,[Base.sum(a[n-i]*b[1+i+j+m-n] for i=0:n-1-j) for j=1:min(n-1,p-m+1)])
        p+1 == length(o) && return o
    else
        o = [Base.sum(b[1+j-i]*a[1+i] for i=0:j) for j=0:min(m-1,p)]
        p+1 == length(o) && return o
        Base.append!(o,[Base.sum(b[m-i]*a[1+i+j] for i=0:m-1) for j=1:min(n-m,p-m+1)])
        p+1 == length(o) && return o
        Base.append!(o,[Base.sum(b[m-i]*a[1+i+j+n-m] for i=0:m-1-j) for j=1:min(m-1,p-n+1)])
        p+1 == length(o) && return o
    end

    return o

end
