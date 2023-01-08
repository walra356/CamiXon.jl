# SPDX-License-Identifier: MIT

# =================================== VectorRational ====================================

@doc raw"""
    VectorRational{T}

Object to decompose a vector of rational numbers

The fields are:
* `.num::Vector{Int}``: vector of normalized numerators
* `.den::Int`: common denominator
* `.val::Vector{Rational}`: vector of rational numbers (simplified = not normalized)
"""
struct VectorRational{T}

    num::Vector{T}
    den::T
    val::Vector{Rational{T}}

end

# ==================================== castVectorRational(vec) =======================

@doc raw"""
    castVectorRational(vec::Vector{Rational{T}}) where T<:Union{Int,BigInt}

Decompose vector of rational numbers.
#### Example:
```
v = [2//3,4//5]
castVectorRational(v)
    VectorRational([10, 12], 15, Rational{Int64}[2//3, 4//5])
```
"""
function castVectorRational(vec::Vector{Rational{T}}) where T<:Union{Int,BigInt}

    val = Base.gcd(vec)
    den = val.den
    num = Base.convert(Vector{T}, (vec .* den))

    return VectorRational(num, den, vec)

end

# ==================================== bernoulli_numbers(nmax) =================

function _bn(o::Vector{Rational{T}}, nstart::Int, nstop::Int) where {T<:Integer}

    for n = nstart:nstop
        a = 0
        if Base.isodd(n)
            b = 1
            for j = 1:n-1
                a -= o[j] * b
                b *= (n + 1 - j)
                b ÷= j        # binomial coefficients are integers
            end
        end
        push!(o, a // n)
    end

    return o

end
# ..............................................................................
function _bn_Int(nmax::Int, nc::Int)

    nmax == 1 && return [1 // 1]

    o = [1 // 1, -1 // 2]
    o = _bn(o, 3, min(nmax, nc))

    return o

end
# ..............................................................................
function _bn_BigInt(o::Vector{T}, nmax::Int, nc::Int) where {T<:Real}

    nmax > nc || return o

    o = convert(Vector{Rational{BigInt}}, o)
    o = _bn(o, nc + 1, nmax)

    return o

end
# ..............................................................................

@doc raw"""
    bernoulli_number(n::T [; msg=true]) where {T<:Integer}
    bernoulli_numbers(nmax::T [; msg=true]) where {T<:Integer}

Bernoulli number array, ``[B_0,⋯\ B_{nmax}]``. This is the *even index convention* 
(odd-indexed numbers vanish - except B[1]=-1/2). The array is calculated by 
repetative use of the recurrence relation
```math
    B_n = - \frac{1}{n+1}\sum_{k=0}^{n-1}\frac{(n+1)!}{k!(n+1-k)}B_k.
```
Special numbers: ``B_0=1,\ B_1=-1/2,\ B_{2n+1}=0\ (\rm{for}\ n>1)``. Starting 
at ``B_0`` is called the *even index convention*.

Integer-overload protection: for `n > 35` (`nmax > 36`) the output is 
autoconverted to `Rational{BigInt}`. *Option*: warning message displayed when 
autoconversion is activated (default: no message).
### Examples:
```
bernoulli_number(60)
#  Warning: output converted to BigInt
#
#  -1215233140483755572040304994079820246041491//56786730

nmax = 10
o = bernoulli_numbers(nmax); println(o)
#  Rational{Int64}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1, 1//42, 0//1, -1//30, 0//1]

n = 60
bernoulli_number(n) == bernoulli_numbers(n+1)[end]             # unit based array
#  true
```
"""
function bernoulli_numbers(nmax::T; msg=true) where {T<:Integer}

    nmax ≥ 0 || return nothing

    nc = 36          # n > nc: integer overload

    V = ConditionalType(nmax, nc; msg)

    nmax = convert(Int, nmax)
    B = _bn_Int(nmax, nc)
    B = _bn_BigInt(B, nmax, nc)

    if (T == BigInt) & (V == Int)
        B = convert(Vector{Rational{BigInt}}, B)
    end

    return B

end
function bernoulli_number(n::T; msg=true) where {T<:Integer}

    n == 0 && return Rational{T}(1 // 1)

    B = bernoulli_numbers(n + 1; msg)[end]

    return B

end

# ==================================== factorialbig(n) =========================

@doc raw"""
    factorialbig(n::Int)

The product of all *positive* integers less than or equal to `n`,
```math
!(n)=n(n-1)(n-2)⋯1.
```
By definition
```math
!(0)=1
```
For *negative* integers the factorial is zero.
#### Examples:
```
factorialbig(20)==factorial(20)
    true

factorialbig(21)
    51090942171709440000

factorial(21)
    OverflowError: 21 is too large to look up in the table; consider using `factorial(big(21))` instead
```
"""
function factorialbig(n::Int)

    n > 20 || return factorial(n)

    return factorial(big(n))

end

# ==================================== faulhaber_polynom(p) ====================

@doc raw"""
    faulhaber_polynom(p [, T=Int])

Vector representation of the Faulhaber polynomial of degree ``p``,
```math
    F(n,p)=\frac{1}{p}\sum_{j=1}^{p}{\binom {p}{p-j}}B_{p-j}n^{j}.
```
``F(n,p)=`` `polynomial(c,n)`, where ``c=[c_0,⋯\ c_p]`` is the coefficient vector, with
```math
    c_0=0,\ c_j=\frac{1}{p}{\binom {p}{p-j}}B_{p-j},
```
with ``j∈\{ 1,⋯\ p\}``. The ``B_0,⋯\ B_{p-1}`` are Bernoulli numbers
(but with ``B_1=+\frac{1}{2}`` rather than ``-\frac{1}{2}``).
### Example:
```
faulhaber_polynom(6)
7-element Vector{Rational{Int64}}:
  0//1
  0//1
 -1//12
  0//1
  5//12
  1//2
  1//6
```
"""
function faulhaber_polynom(k::Int; T=Int)

    k < 1 && return 0
    k > 1 || return 1//1

    P = CamiXon.pascal_triangle(k)[end][1:end-1]
    B = CamiXon.bernoulli_numbers(k); B[2]=-B[2]  # was bernoulli_numbers(k-1; T)

    F = (B .* P)  // k

    F = Base.append!(F,0//1)   # add polynomial constant (zero in this case)

    return Base.reverse(F)     # reverse to standard order

end
function faulhaber_polynom1(k::Int)       # short argument: better performance

    k < 1 && return 0
    k > 1 || return 1//1

    P = CamiXon.pascal_triangle(k)[end][1:end-1]
    B = CamiXon.bernoulli_numbers(k-1); B[2]=-B[2]

    F = (B .* P)  // k

    F = Base.append!(F,0//1)   # add polynomial constant (zero in this case)

    return Base.reverse(F)     # reverse to standard order

end

# =================================== faulhaber_summation(n,p;T) ===============

@doc raw"""
    faulhaber_summation(n, p [, T=Int])

Sum of powers of natural numbers ``1,⋯\ n``,
```math
    FS(n,p)=\sum_{k=1}^{n}k^{p}=F(n,p+1).
```
where ``F(n,p)`` is the Faulhamer polynomial of degree ``p``.
### Examples:
```
faulhaber_summation(5,1)
 15

faulhaber_summation(3,60; T=BigInt)
  42391158276369125018901280178
```
"""
function faulhaber_summation(n::Int, p::Int; T=Int)

    n ≠ 0 || return nothing

    F = CamiXon.faulhaber_polynom(p+1; T)
    o = 0
    for k=1:p+1
        for i=1:k
            F[k+1] *= n # avoid n^k in o = Base.sum([F[k+1]*n^k for k=1:p+1])
        end
        o += F[k+1]
    end

    Base.denominator(o) == 1 || error("Error: Faulhaber sum failed")

    return Base.numerator(o)

end
function faulhaber_summation(n::Int, p::Int)   # short argument: better performance

    n ≠ 0 || return 0

    F = CamiXon.faulhaber_polynom(p+1)
    o = 0
    for k=1:p+1
        for i=1:k
            F[k+1] *= n # avoid n^k in o = Base.sum([F[k+1]*n^k for k=1:p+1])
        end
        o += F[k+1]
    end

    Base.denominator(o) == 1 || error("Error: Faulhaber sum failed")

    return Base.numerator(o)

end

# ============================= Fibonacci numbers =================================



# ======================== _collect_Int and _collect_BigInt ====================

function _fn(o, nstart::Int, nstop::Int)

    for n = nstart:nstop
        push!(o, o[n-1] + o[n-2])
    end

    return o

end
# ..............................................................................
function _fn_Int(o, nstop::T, nc::Int) where {T<:Integer}

    nstop = convert(Int, nstop)

    nstart = length(o) + 1
    nstop = min(nstop, nc)

    o = _fn(o, nstart, nstop)

    return o

end
# ..............................................................................
function _fn_BigInt(o, nstop::T, nc::Int) where {T<:Real}

    nstop = convert(Int, nstop)

    nstop > nc || return o

    o = bigconvert(o)
    o = _fn(o, nc + 1, nstop)

    return o

end
# ..............................................................................

@doc raw"""
    fibonacci_number(n::T [; msg=true]) where T<:Integer
    fibonacci_numbers(nmax::T [; msg=true]) where T<:Integer

A sequence of integers,  ``F_0,⋯\ F_{nmax}``, in which each element is the sum of the 
two preceding ones, 
```math
    F_n = F_{n-1}+F_{n-2}.
```
where ``F_1=1`` and ``F_0=0``. 

Integer-overload protection: for `nmax > 92` the output is autoconverted to BigInt.
Optional: a warning message is displayed when autoconversion is activated (default: no message)
#### Example:
```
julia> fibonacci_numbers(20)
21-element Vector{Int64}:
   0
   1
   1
   2
   3
   5
   8
  13
   ⋮
1597
2584
4181

julia> fibonacci_number(100)
Warning: converted to BigInt (integer overload protection)
354224848179261915075

julia>
```
"""
function fibonacci_numbers(nmax::T; msg=true) where {T<:Integer}

    nmax = nmax + 1

    nmax > 0 || return nothing
    nmax > 1 || return T(0)

    nc = 93   # n > nc: integer overload protection activated

    o = _fn_Int([0, 1], nmax, nc)
    o = _fn_BigInt(o, nmax, nc)
    o = ((T == BigInt) & !protectInt(nmax, nc; msg)) ? bigconvert(o) : o

    return o

end
# ..............................................................................
function fibonacci_number(n::T; msg=true) where {T<:Integer}

    n > 0 || return T(0)

    F = fibonacci_numbers(n; msg)[end]

    return F

end

# =================================== harmonic number(n;T) ===============

function _hn(o, nstart::Int, nstop::Int)

    for n = nstart:nstop
        o += 1 // n
    end

    return o

end
# ..............................................................................
function _hn_Int(o, nstop::T, nc::Int) where {T<:Integer}

    nstop = convert(Int, nstop)

    nstart = length(o) + 1
    nstop = min(nstop, nc)

    o = _fn(o, nstart, nstop)

    return o

end
# ..............................................................................
function _hn_BigInt(o, nstop::T, nc::Int) where {T<:Real}

    nstop = convert(Int, nstop)

    nstop > nc || return o

    o = bigconvert(o)
    o = _hn(o, nc + 1, nstop)

    return o

end
# ..............................................................................

@doc raw"""
    harmonic_number(n::T [; msg=true]) where {T<:Integer} 

Sum of the reciprocals of the first ``n`` natural numbers
```math
    H_n=\sum_{k=1}^{n}\frac{1}{k}.
```
Integer-overload protection: for `n > 46` the output is autoconverted to Rational{BigInt}.
Optional: a warning message is displayed when autoconversion is activated (default: no message)
### Examples:
```
o = [harmonic_number(i) for i=1:7]; println(o)
#  Rational{Int64}[1//1, 3//2, 11//6, 25//12, 137//60, 49//20, 363//140]

o = [harmonic_number(46; msg=true)]; println(o)
#  Rational{Int64}[5943339269060627227//1345655451257488800]

o = [harmonic_number(47; msg=true)]; println(o)
#  Warning: output converted to BigInt
#  Rational{BigInt}[280682601097106968469//63245806209101973600]

o = [harmonic_number(47)]; println(o)
#  Rational{BigInt}[280682601097106968469//63245806209101973600]

harmonic_number(60; msg=false)
#  15117092380124150817026911//3230237388259077233637600

harmonic_number(12) == harmonic_number(12, 1)
#  true
```
"""
function harmonic_number(n::T; msg=true) where {T<:Integer}

    n ≠ 0 || return 0

    nc = 46

    o = _fn_Int(1//1, nmax, nc)
    o = _fn_BigInt(o, nmax, nc)
    o = (T == BigInt) & !protectInt(nmax, nc; msg) ? bigconvert(o) : o

    return o

end

# =================================== harmonic number(n, p [; msg=false]) ===============

function _hn_Int(n, nc, p)

    n = Int(n)

    o = 0 // 1
    for j = 1:min(n, nc)
        a = 1
        for i = 1:p
            a *= j
        end
        o += 1 // a
    end

    return o

end
# .......................................................................................
function _hn_BigInt(o, n, nc, p)

    one = big(1)

    for j = nc+1:n
        a = one
        for i = 1:p
            a *= big(j)
        end
        o += one // a
    end

    return o

end
@doc raw"""
    harmonic_number(n::T, p::Int [; msg=false]) where {T<:Integer}

Sum of the ``p_{th}`` power of reciprocals of the first ``n`` numbers
```math
    H_{n,p}=\sum_{k=1}^{n}\frac{1}{k^p}.
```
Integer-overload protection: the output is autoconverted to Rational{BigInt} when appropriate.
Optional: a warning message is displayed when autoconversion is activated (default: no message)
### Examples:
```
o = [harmonic_number6(46,1; msg=true)]; println(o)
#  Rational{Int64}[5943339269060627227//1345655451257488800]

o = [harmonic_number6(47,1; msg=true)]; println(o)
# Warning: output converted to BigInt
# Rational{BigInt}[280682601097106968469//63245806209101973600]

o = [harmonic_number6(47,1)]; println(o)
# Rational{BigInt}[280682601097106968469//63245806209101973600]

harmonic_number(12, -3) == faulhaber_summation(12, 3)
  true
```
"""
function harmonic_number(n::T, p::Int; msg=false) where {T<:Integer}

    n ≠ 0 || return 0

    ac = [46, 24, 16, 12, 10, 8, 7, 6, 6, 6, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3]
    nc = p > 17 ? 0 : ac[p]

    if p > 0
        o = _hn_Int(n, nc, p)
        if n > nc
            T ≠ BigInt ? (msg && println("Warning: output converted to BigInt")) : false
            o *= big(1)
            o = _hn_BigInt(o, n, nc, p)
        else
            T ≠ BigInt ? true : o *= T(1)
        end
    else
        p = -p
        F = CamiXon.faulhaber_polynom(p + 1; T)
        o = 0
        for k = 1:p+1
            for i = 1:k
                F[k+1] *= n
            end
            o += F[k+1]
        end
        Base.denominator(o) == 1 || error("Error: Faulhaber sum failed")
        o = Base.numerator(o)
    end

    return o

end
function harmonic_number1(n::Int, p::Int)    # short argument: better performance

    n ≠ 0 || return 0

    if p > 0
        o = 0//1
        for j=1:n
            a = 1
            for i=1:p
                a *= j
            end
        o += 1//a
        end
    else
        p = -p
        F = CamiXon.faulhaber_polynom(p+1)
        o = 0
        for k=1:p+1
            for i=1:k
                F[k+1] *= n
            end
            o += F[k+1]
        end
        Base.denominator(o) == 1 || error("Error: Faulhaber sum failed")
        o = Base.numerator(o)
    end

    return o

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

# ==================================== pascal_triangle(nmax)  ============

@doc raw"""
    pascal_triangle(nmax [, T=Int])

Pascal triangle of binomial coefficients ``\binom{n}{k}`` for ``n=0,\ 1,⋯\ nmax``
### Example:
```
pascal_triangle(5)
6-element Vector{Vector{Int64}}:
 [1]
 [1, 1]
 [1, 2, 1]
 [1, 3, 3, 1]
 [1, 4, 6, 4, 1]
 [1, 5, 10, 10, 5, 1]
```
"""
function pascal_triangle(nmax::T) where {T<:Integer}

    nmax < 0 && error("Error: nmax must be a non-negative integer")
    nmax > T(10000) && error("Error: integer overflow")

    o = [Base.ones(T, n + 1) for n = 0:nmax]

    for n = 2:nmax
        for k = 1:n÷2
            o[n+1][k+1] = o[n][k+1] + o[n][k]
            o[n+1][n+1-k] = o[n+1][k+1]
        end
    end

    return o

end

# ==================================== pascal_next(nmax)  ======================

@doc raw"""
    pascal_next(nmax)

Next row of Pascal triangle
### Example:
```
a = [1, 4, 6, 4, 1]
pascal_next(a)
 [1, 5, 10, 10, 5, 1]
```
"""
function pascal_next(a::Vector{Int})

    n = Base.length(a) + 1
    o = Base.ones(Int,n)

    for k=1:n÷2
        o[k+1] = a[k+1] + a[k]
        o[n-k] = o[k+1]
    end

    return o

end

# ====================== permutations_unique_count(p, i) =======================

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

# ===================== Pochhammer(x, p) =======================================

@doc raw"""
    pochhammer(x::T, p::Int) where T<:Real

Pochhammer symbol ``(x)_{p}`` for integral ``p``,
```math
(x)_{p}=\begin{cases}
1 & p=0\\
x(x+1)(x+2)⋯(x+p-1) & p>0
\end{cases}
```

Note that ``(x)_{p}=0`` for ``x=0,-1,⋯\ -(p-1)``
#### Examples:
```
x = [-4,-3,-2,-1, 0, 1, 2 , 3, 4]
pochhammer.(x,5) == [0, 0, 0, 0, 0, 120, 720, 2520, 6720]
  true

pochhammer.(x,0) == [1, 1, 1, 1, 1, 1, 1, 1, 1]
  true

o = [pochhammer.([x for x=0:-1:-p],p) for p=0:5]
println("non-positive integer x = 0,⋯\ -p:")
for p=0:5
    println("p = $p: $(o[p+1])")
end
  non-positive integer x = 0,⋯\ -p:
  p = 0: [1]
  p = 1: [0, -1]
  p = 2: [0, 0, 2]
  p = 3: [0, 0, 0, -6]
  p = 4: [0, 0, 0, 0, 24]
  p = 5: [0, 0, 0, 0, 0, -120]

 o = [pochhammer.([x for x=0:p],p) for p=0:5]
 println("non-negative integer x = 0,⋯\  p:")
 for p=0:5
     println("p = $p: $(o[p+1])")
 end
   non-negative integer x = 0,⋯\  p:
   p = 0: [1]
   p = 1: [0, 1]
   p = 2: [0, 2, 6]
   p = 3: [0, 6, 24, 60]
   p = 4: [0, 24, 120, 360, 840]
   p = 5: [0, 120, 720, 2520, 6720, 15120]

x = -1//50
pochhammer(x,20)
  OverflowError: -1491212300990613201 * 449 overflowed for type Int64

x = convert(Rational{BigInt}, -1//50)
pochhammer(x,20)
  -21605762356630090481082546653745369902321614221999//9536743164062500000000000000000000
```
"""
function pochhammer(x::T, p::Int) where T<:Real

    p > 0 || return 1

    o = x

    for n=1:p-1
        o *= (x+n)
    end

    return o

end
# ============================== triangle_coefficient(a, b, c) =============================

@doc raw"""
    triangle_coefficient(a::Real, b::Real, c::Real)

Triangle coefficient for a triangle of sides `a`, `b` and `c`.

#### Example:
```
triangle_coefficient(3, 4, 5)
    1//180180

triangle_coefficient(1//2, 1, 1.5)
    1//12
```
"""
function triangle_coefficient(a::Real, b::Real, c::Real)

    (a,b,c) = promote(a,b,c)

    isinteger(a + b + c) || return 0

    A = Int(a + b - c)
    B = Int(b + c - a)
    C = Int(c + a - b)

    A = A ≥ 0 ? factorialbig(A) : return 0
    B = B ≥ 0 ? factorialbig(B) : return 0
    C = C ≥ 0 ? factorialbig(C) : return 0

    num = A * B * C
    den = factorialbig(Int(a+b+c+1))

    return num//den

end

# ============================ istriangle(a, b, c) =============================

@doc raw"""
    istriangle(a::Real, b::Real, c::Real)

Triangle condition for a triangle of sides `a`, `b` and `c`.

#### Example:
```
istriangle(3, 4, 5)
    true

istriangle(1//2, 1, 1.5)
    true
```
"""
function istriangle(a::Real, b::Real, c::Real)

    Δ = triangle_coefficient(a,b,c)

    valid = Δ > 0 ? true : false

    return valid

end

# ...................... texp(x, p) .........................................

function _texp_int(x, p::Int)

    o = y = typeof(x)(1)

    x ≠ 0 || return o

    for n=1:p
        y *= x//n
        o += y
    end

    return o

end

function _texp_real(x, p::Int)

    o = y = typeof(x)(1)

    x ≠ 0.0 || return o

    for n=1:p
        y *= x/n
        o += y
    end

    return o

end

@doc raw"""
    texp(x::T, a::T, p::Int) where T <: Real

Taylor expansion of exp(x) about ``x = a`` up to order p.
```math
    \mathsf{texp}(x,a,p) = 1 + (x-a) + \frac{1}{2}(x-a)^2 + ⋯ + \frac{1}{p!}(x-a)^p.
```
### Examples:
```
p = 5
texp(1.0, 0.0, 5)
 2.7166666666666663

texp(1, 0, 5)
 163//60
```
"""
function texp(x::T, a::T, p::Int) where T <: Real

    x = x - a

    V = typeof(x)

    return  V <: Rational ? _texp_int(x, p) : V <: Integer ? _texp_int(x, p) : _texp_real(x, p)

end
