

# ...................................................... VectorRational .........................................................

@doc raw"""
    VectorRational

Object to decompose a vector of rational numbers

The fields are:
* `.num::Vector{Int}``: vector of normalized numerators
* `.den::Int`: common denominator
* `.val::Vector{Rational}`: vector of rational numbers (simplified = not normalized)
"""
struct VectorRational

    num::Vector{Int}
    den::Int
    val::Vector{Rational{Int64}}

end

# ==================================== analyzeVectorRational(vec) =======================

@doc raw"""
    normalize_VectorRational(vec::Vector{Rational{Int}})

Decompose vector of rational numbers.
#### Example:
```
v = [2//3,4//5]
normalize_VectorRational(v)
 VectorRational([10, 12], 15, Rational{Int64}[2//3, 4//5])
```
"""
function normalize_VectorRational(vec::Vector{Rational{Int}})

    val = gcd(vec)
    den = val.den
    num = convert(Vector{Int},(vec .* den))

    return VectorRational(num, den, vec)

end

# ==================================== bernoulli_numbers(nmax) =================

@doc raw"""
    bernoulli_numbers(nmax)

Bernoulli numbers ``B_0,\ \cdots,\ B_{nmax}`` calculated by repetative use of the recurrence relation
```math
    B_n = - \frac{1}{n+1}\sum_{k=0}^{n-1}\frac{(n+1)!}{k!(n+1-k)}B_k
```
Special numbers: ``B_0=1,\ B_1=-1/2,\ B_{2n+1}=0\ (\rm{for}\ n>1)``.
### Examples:
```
bernoulli_numbers(10)
11-element Vector{Rational{Int64}}:
  1//1
 -1//2
  1//6
  0//1
 -1//30
  0//1
  1//42
  0//1
 -1//30
  0//1
  5//66
```
"""
function bernoulli_numbers(nmax::Int)

    B = nmax < 36 ? zeros(Rational{Int},nmax+1) : zeros(Rational{BigInt},nmax+1)

    B[1] = 1//1
    B[2] = -1//2

    for m = 3:nmax+1
        B[m] = 0//1
        if isodd(m)
            b = 1//m
            for j = 1:m-1
                B[m] -= b * B[j]
                b *= (m+1-j) // j      # binomial coefficients
            end
        end
    end

    return B

end

# ==================================== _canonical_partition(n, m) =======================

function _canonical_partition(n::Int, m::Int)

    o = Base.fill(m,Base.cld(n,m))                              # init partition
    o[Base.cld(n,m)]=((n%m)≠0 ? n%m : m)                        # adjust last element of partition

    return o

end

"""
    canonical_partitions(n; header=false, reverse=true)

The canonical partition in integers of the integer n

header=true : unit patition included in output
#### Examples:
```
canonical_partitions(6; header=true, reverse=false)
6-element Array{Array{Int64,1},1}:
 [6]
 [5, 1]
 [4, 2]
 [3, 3]
 [2, 2, 2]
 [1, 1, 1, 1, 1, 1]

canonical_partitions(6; header=true)
6-element Array{Array{Int64,1},1}:
 [1, 1, 1, 1, 1, 1]
 [2, 2, 2]
 [3, 3]
 [4, 2]
 [5, 1]
 [6]

canonical_partitions(6)
5-element Array{Array{Int64,1},1}:
 [1, 1, 1, 1, 1, 1]
 [2, 2, 2]
 [3, 3]
 [4, 2]
 [5, 1]
```
"""
function canonical_partitions(n::Int, m=0; header=true, reverse=true)

    h = header ? n : n-1

    if m == 0
        if reverse
            o = [_canonical_partition(n,m) for m=1:h]
        else
            o = [_canonical_partition(n,m) for m=h:-1:1]
        end
    elseif 0 < m <= n
        o = _canonical_partition(n,m)
    else
        o = nothing
    end

    return o

end



function _partition_count(n::Int,k::Int)

    (n<0)|(k<0)|(k>n) ? 0 : (k==n)|(k==1) ? 1 : _partition_count(n-k,k) + _partition_count(n-1,k-1)

end

function _partition(a::Array{Int,1}, n::Int, i::Int, cp::Array{Array{Array{Int,1},1},1})

    o = a[1:i-1]
    m = a[i]-1                                           # m: partition value
    ni = n - Base.sum(o)                                 # ni: sub-partition index at partition index i

    Base.append!(o,cp[ni][m])                            # complete partition by appending it to a

    return o

end

function _restricted_partitions(o::Array{Int,1}, n::Int, np::Int, cp::Array{Array{Array{Int,1},1},1})

    oo = [o]

    for p=1:np-1
        i = Base.findlast(x -> x > 1, oo[p])
        Base.append!(oo,[_partition(oo[p],n,i,cp)])
    end

    return oo

end

"""
    integer_partitions(n [,m]; transpose=false, count=false)

default              : The integer partitions of n

count=true           : The number of integer partitions of n

transpose=false/true : for m>0 restricted to partitions with maximum part/length m

definitions:

The integer partition of the positive integer n is a nonincreasing sequence of positive integers p1, p2,... pk whose sum is n.

The elements of the sequence are called the parts of the partition.
#### Examples:
```
integer_partitions(7)
15-element Array{Array{Int64,1},1}:
 [1, 1, 1, 1, 1, 1, 1]
 [2, 2, 2, 1]
 [3, 3, 1]
 [4, 3]
 [5, 2]
 [6, 1]
 [7]
 [2, 2, 1, 1, 1]
 [3, 2, 2]
 [4, 2, 1]
 [5, 1, 1]
 [2, 1, 1, 1, 1, 1]
 [3, 2, 1, 1]
 [4, 1, 1, 1]
 [3, 1, 1, 1, 1]

integer_partitions(7; count=true)
15

integer_partitions(7,4; count=true)
3

integer_partitions(7,4)
3-element Array{Array{Int64,1},1}:
 [4, 3]
 [4, 2, 1]
 [4, 1, 1, 1]

integer_partitions(7,4; transpose=true)
3-element Array{Array{Int64,1},1}:
 [2, 2, 2, 1]
 [3, 2, 1, 1]
 [4, 1, 1, 1]
```
"""
function integer_partitions(n::Int, m=0; transpose=false, count=false)

    cp = [canonical_partitions(m) for m=1:n]
    pc = [_partition_count(n,m)  for m=1:n]
    oo = [ones(Int,n)]

    np = m > 0 ? pc[m] : sum(pc)

    if !count

        if m == 0
            o = [_restricted_partitions(cp[n][p],n,pc[p],cp) for p=2:n]
            for p=1:n-1 append!(oo,o[p]) end
        else
            oo = _restricted_partitions(cp[n][m],n,pc[m],cp)
        end

        if transpose
            for p=1:np
                l = length(oo[p])
                s=max(oo[p][1],l)
                mat = zeros(Int,s,s)
                for j=1:l for i=1:oo[p][j] mat[i,j]=1 end end
                oo[p] = [sum(mat[i,:]) for i=1:oo[p][1]]
            end

        end

    end

    return count ? np : oo

end

# ===================================== log10_characteristic_power(x) ==============================================
"""
    log10_characteristic_power(x)

characteristic power-of-10 of the number x
#### Examples:
```
log10_characteristic_power.([3,30,300])
3-element Vector{Int64}:
 0
 1
 2
```
"""
log10_characteristic_power(x) = Base.round(Int,Base.floor(log10(x)))

# ==================================== log10_mantissa(x) ============================================================
"""
    log10_mantissa(x)

log10 mantissa of the number x
#### Examples:
```
log10_mantissa.([3,30,300])
3-element Vector{Float64}:
 0.47712125471966244
 0.4771212547196624
 0.4771212547196626
```
"""
log10_mantissa(x) = Base.log10(x)-Base.floor(Base.log10(x))

# ==================================== polynom(c,x) ============================================================

@doc raw"""
    polynom(coords,x)

Method to evaluate the function ``f(x)=polynom(c,x)``, where
``c=[c_0,\ \ldots,\ c_d]`` is the vector representation of a polynomial of degree ``d``.
```math
    p(c,x)=c_0 + c_1 x + \cdots + c_d x^d.
```
### Examples:
```
coords = ones(6)             # for polynomial of degree 5 with unit coefficients
f(x) = polynom(coords,x)
println([f(1.0),f(2.0)])     # values of polynomial for x = 1.0 and x = 2.0
 [6.0, 63.0]
```
"""
function polynom(coords::Vector{T}, x::T) where T<:Number

    k = Base.length(coords)
    X = Base.ones(T,k)

    for i=2:k
        X[i] = X[i-1] * x
    end

    return LinearAlgebra.dot(coords, X)

end

# ==================================== polynom_derivative(coords) ==============

@doc raw"""
    polynom_derivative(coords)

Vector representation of the first derivative of the polynomial `coords`,
```math
    p'(c,x)=c_1 + 2 c_2 x + \cdots + d c_d x^{d-1},
```
Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``.
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,\ \ldots,\ c_d]``
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
function polynom_derivative(coords::Vector{<:Number})

    k = Base.length(coords)
    k > 1 || return [0]

    return coords[2:end] .* Base.OneTo(k-1)

end

# =============================== polynom_derivative(coords[,deriv=0]) =========

@doc raw"""
    polynom_derivatives(coords[,deriv=0])

Vector representation of derivatives of the polynomial `coords`.

Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``.
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,\ \ldots,\ c_d]``
consisting of the polynomial coefficients.

`deriv`: derivative of choice; `default`: collection of all (nontrivial) derivatives.

### Examples:
```
coords=[1,1,1,1,1]               # vector representation of a polynomial of degree d=4
polynom_derivatives(coords)      # `all' (nontrivial) derivatives of polynomial `coords`
5-element Vector{Vector{Int64}}:
 [1, 2, 3, 4]
 [2, 6, 12]
 [6, 24]
 [24]
 [0]

polynom_derivatives(coords; deriv=2)          # second derivative of polynomial `coords`
3-element Vector{Int64}:
  2
  6
 12
```
"""
function polynom_derivatives(coords::Vector{<:Number}; deriv=0)

    deriv < 0 && error("jwError: negative derivative not defined")

    k = deriv > 0 ? deriv+1 : Base.length(coords)

    coords = CamiXon.polynom_derivative(coords)

    deriv ≠ 1 ? o = [coords] : return coords

    for i=2:k-1
        coords = polynom_derivative(coords)
        Base.push!(o,coords)
    end

    deriv == 0 ? Base.push!(o,[0]) : return o[deriv]

    return o

end

# ==================================== polynom_power(coords, p) ================

@doc raw"""
    polynom_power(coords, p)

Vector representation of the polynomial `coords` raised to the power `p` which
results in a polynomial in a vector space of dimension ``p d + 1``.

Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``.
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,\ \ldots,\ c_d]``
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
function polynom_power(coords::Vector{<:Number}, power::Int)

    power >= 0 || error("jwError: negative powers not allowed")
    power == 2 && return polynom_product(coords, coords)
    power == 1 && return coords
    power == 0 && return [1]

    o = polynom_product(coords, coords)

    for i=1:power-2
        o = polynom_product(o, coords)
    end

    return o

end

# ==================================== polynom_powers(coords, pmax) ============

@doc raw"""
    polynom_powers(coords, pmax)

The polynomial `coords` raised to the powers 1,...,pmax  which
results in a collection of polynomials in vector spaces of dimension ``d+1`` tot ``p d + 1``.

Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``.
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,\ \ldots,\ c_d]``
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
function polynom_powers(coords::Vector{<:Number}, pmax::Int)

    pmax > 0 || error("jwError: minimum power included is unity")

    o = [coords]

    for i=1:pmax-1
        push!(o,polynom_product(o[end], coords))
    end

    return o

end

# ==================================== polynom_primitive(coords) ====================

@doc raw"""
    polynom_primitive(coords)

Vector representation of the primitive of the polynomial `coords` which is a
polynomial in a vector space of dimension ``p d + 1``.
```math
    P(c,x)=c_{int} +c_0 x + \frac{1}{2} c_1 x^2 + \frac{1}{3} c_2 x^3 + \cdots + \frac{1}{d+1} c_d x^{d+1},
```
The constant of integration is set to zero, ``c_{int} = 0``.

Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``.
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,\ \ldots,\ c_d]``
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
function polynom_primitive(coords::Vector{<:Number})

    d = [1//p for p ∈ Base.eachindex(coords)]

    coords = coords .* d

    return Base.pushfirst!(coords,0)      # constant of integration equal to zero

end

# ==================================== polynom_product(a, b) ============================================================

@doc raw"""
    polynom_product(a::Vector{<:Number}, b::Vector{<:Number})

Vector representation of the product of two polynomials, ``a`` and ``b`` which
is a polynomial in a vector space of dimension ``d=m+n``,
```math
    p(c,x)=a_0b_0 + (a_0b_1 + b_0a_1)x + \cdots + a_n b_m x^{n+m}.
```
Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``.
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,\ \ldots,\ c_d]``
consisting of the polynomial coefficients.
####
```
a = [1,1]
b = [1,- 1]
o = polynomial_product(a, b); println(o)
 [1, 0, -1]
```
"""
function polynom_product(a::Vector{<:Number}, b::Vector{<:Number})

    n = Base.length(a)
    m = Base.length(b)

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
    polynom_product_expansion(a::Vector{<:Number}, b::Vector{<:Number}, p::Int)

Vector representation of the product of two polynomials, ``a`` (degree=``n``) and ``b`` (degree=``m``), with ``m≤n``
truncated at degee=``m`` is a polynomial in a vector space of dimension ``d=m+1``. If ``p`` is the `polynom_product,
the polynom_product_expansion is ``p[1:p+1]``
####
```
a = [1,-1,1]
b = [1,1,-1]
o = polynom_product(a, b); println(o)
 [1, 0, -1, 2, -1]

o = polynom_product_expansion(a, b, 3); println(o)
 [1, 0, -1]
```
"""
function polynom_product_expansion(a::Vector{<:Number}, b::Vector{<:Number}, p::Int)

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

# ==================================== permutations_unique_count(p, i) =======================

@doc raw"""
    permutations_unique_count(p::Array{Array{Int64,1},1}, i::Int)

Number of unique permutations of the subarray ``p[i]``.
#### Example:
```
p = [[1,2,3],[2,3,1,4,3]]
permutations_unique_count(p,2)
 60
```
"""
function permutations_unique_count(p::Array{Array{Int64,1},1}, i::Int)

    o = Base.factorial(Base.length(p[i]))
    d = Base.Dict([(n,Base.count(x->x==n,p[i])) for n ∈ Base.unique(p[i])])

    for j ∈ Base.eachindex(Base.unique(p[i]))
        o = o ÷ Base.factorial(d[Base.unique(p[i])[j]])
    end

    return o

end
