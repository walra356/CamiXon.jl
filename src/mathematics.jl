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
julia> v = [2//3,4//5];
julia> castVectorRational(v)
VectorRational([10, 12], 15, Rational{Int64}[2//3, 4//5])
```
"""
function castVectorRational(vec::Vector{Rational{T}}) where T<:Union{Int,BigInt}

    val = Base.gcd(vec)
    den = val.den
    num = Base.convert(Vector{T}, (vec .* den))

    return VectorRational(num, den, vec)

end

# ============================ bernoulliB_array ==============================

global glBn_Int = [1 // 1, -1 // 2, 1 // 6, 0 // 1, -1 // 30, 0 // 1, 1 // 42,
    0 // 1, -1 // 30, 0 // 1, 5 // 66, 0 // 1, -691 // 2730, 0 // 1,
    7 // 6, 0 // 1, -3617 // 510, 0 // 1, 43867 // 798, 0 // 1,
    -174611 // 330, 0 // 1, 854513 // 138, 0 // 1, -236364091 // 2730,
    0 // 1, 8553103 // 6, 0 // 1, -23749461029 // 870, 0 // 1,
    8615841276005 // 14322, 0 // 1, -7709321041217 // 510, 0 // 1,
    2577687858367 // 6, 0 // 1]

global glBn_BigInt = convert(Vector{Rational{BigInt}}, glBn_Int)


# ..............................................................................
#function _bn_Int(n::Int, nc::Int)             kanweg ?
#
#    nstop = min(n, nc)

#    o = glBn_Int[1:1+nstop]

#    return o

#end
# ..............................................................................
function _bn_BigInt(n::Int, nc::Int)

    nul = big(0)
    one = big(1)

    o = glBn_BigInt[1:1+nc]
    for m = nc+2:n+1
        a = nul
        if Base.isodd(m)
            b = one
            for j = 1:m-1
                a -= o[j] * b
                b *= (m + 1 - j)
                b ÷= j                     # binomial coefficients are integers
            end
        end
        Base.push!(o, a // big(m))
    end

    return o

end
# ..............................................................................
@doc raw"""
    bernoulliB(n::T [; msg=true]) where {T<:Integer}
    bernoulliB_array(nmax::T [; msg=true]) where {T<:Integer}

Bernoulli number array, ``[B_0,⋯\ B_{nmax}]``. This is the *even index convention* 
(odd-indexed numbers vanish - except B[1]=-1/2). The array is calculated by 
repetative use of the recurrence relation.
```math
    B_n = - \frac{1}{n+1}\sum_{k=0}^{n-1}\frac{(n+1)!}{k!(n+1-k)}B_k.
```
Special numbers: ``B_0=1,\ B_1=-1/2,\ B_{2n+1}=0\ (\rm{for}\ n>1)``. Starting 
at ``B_0`` is called the *even index convention*.

Integer-overflow protection: for `n > 35` the output is autoconverted to `Rational{BigInt}`. 
By default the capture message is activated: 
"Warning: bernoulliB autoconverted to Rational{BigInt}". 
### Examples:
```
julia> bernoulliB(60)
"Warning: bernoulliB autoconverted to Rational{BigInt}"
-1215233140483755572040304994079820246041491//56786730

julia> bernoulliB_array(10; println(o)
Rational{Int64}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1, 1//42, 0//1, -1//30, 0//1]

julia> n = 60;
julia> bernoulliB(n) == bernoulliB_array(n)[end]             
true
```
"""
function bernoulliB(n::T; msg=true) where {T<:Integer}

    o = CamiXon.bernoulliB_array(n; msg)[end]

    return o

end
# ..............................................................................
function bernoulliB_array(nmax::T; msg=true) where {T<:Integer}

    n = Int(nmax)
    nc = 35

    if n ≤ nc
        o = T == Int ? glBn_Int[1:1+n] : glBn_BigInt[1:1+n]
    else
        o = _bn_BigInt(n, nc)
        msg && T == Int && println("Integer-overflow protection: bernoulliB converted to Rational{BigInt}")
    end

    return o

end

#function bernoulliB1(n::T; msg=true) where {T<:Integer}  #kanweg
#
#    n = Int(n)
#    nc = 35
#
#    if n ≤ nc
#       o = T == Int ? glBn_Int[1:1+n][end] : glBn_BigInt[1:1+n][end]
#    else
#        o = _bn_BigInt(n, nc)[end]
#        msg && T == Int && println("Integer-overflow protection: bernoulliB converted to Rational{BigInt}")
#    end
#
#    return o
#
#end


# ==================================== bigfactorial(n) =========================

@doc raw"""
    bigfactorial(n::Int [; msg=true])

The product of all *positive* integers less than or equal to `n`,
```math
!(n)=n(n-1)(n-2)⋯1.
```
By definition
```math
!(0)=1
```
For *negative* integers the factorial is zero.
Integer-overflow protection: for `n > 20` the output is autoconverted to `BigInt`.
By default the capture message is activated: "Warning: bigfactorial autoconverted to BigInt". 
#### Examples:
```
julia> bigfactorial(20) == factorial(20)
true

julia> bigfactorial(21)
Warning: bigfactorial autoconverted to BigInt
51090942171709440000

julia> bigfactorial(21; msg=false)
51090942171709440000

factorial(21)
    OverflowError: 21 is too large to look up in the table; consider using `factorial(big(21))` instead
```
"""
function bigfactorial(n::T; msg=true) where {T<:Integer}
    
    T == Int || return factorial(n)
    
    n = Int(n)
    nc = 20
    
    o = n > nc ? factorial(big(n)) : o = factorial(n)
    
    msg && n > nc && println("Warning: bigfactorial autoconverted to BigInt")

    return o

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
    B = CamiXon.bernoulliB_array(k-1); B[2]=-B[2]  # was bernoulliB_array(k-1; T)

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

global glFn_Int = [
    
    0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 
    4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811, 514229, 
    832040, 1346269, 2178309, 3524578, 5702887, 9227465, 14930352, 24157817, 
    39088169, 63245986, 102334155, 165580141, 267914296, 433494437, 701408733, 
    1134903170, 1836311903, 2971215073, 4807526976, 7778742049, 12586269025, 
    20365011074, 32951280099, 53316291173, 86267571272, 139583862445, 
    225851433717, 365435296162, 591286729879, 956722026041, 1548008755920, 
    2504730781961, 4052739537881, 6557470319842, 10610209857723, 17167680177565, 
    27777890035288, 44945570212853, 72723460248141, 117669030460994, 
    190392490709135, 308061521170129, 498454011879264, 806515533049393, 
    1304969544928657, 2111485077978050, 3416454622906707, 5527939700884757, 
    8944394323791464, 14472334024676221, 23416728348467685, 37889062373143906, 
    61305790721611591, 99194853094755497, 160500643816367088, 259695496911122585, 
    420196140727489673, 679891637638612258, 1100087778366101931, 1779979416004714189, 
    2880067194370816120, 4660046610375530309, 7540113804746346429

]
# ..............................................................................
global glFn_BigInt = convert(Vector{BigInt}, glFn_Int)

# ..............................................................................
function _fn_BigInt(n::Int, nc::Int)

    o = glFn_BigInt[1:1+nc]
    for k = nc+2:n+1
        Base.push!(o, o[k-1] + o[k-2])
    end

    return o

end
# ..............................................................................

@doc raw"""
    fibonacciF(n::T [; msg=true]) where T<:Integer
    fibonacciF_array(nmax::T [; msg=true]) where T<:Integer

The sequence of integers,  ``F_0,⋯\ F_{nmax}``, in which each element is 
the sum of the two preceding ones, 
```math
    F_n = F_{n-1}+F_{n-2}.
```
with ``F_1=1`` and ``F_0=0``. 

Integer-overflow protection: for `n > 92` the output is autoconverted to BigInt. 
By default the capture message is activated: "Warning: fibonacciF autoconverted to BigInt". 
#### Example:
```
julia> fibonacciF_array(20)
21-element Vector{Int64}:
   0
   1
   1
   2
   3
   5
   ⋮
1597
2584
4181

julia> fibonacciF(100)
Warning: fibonacciF autoconverted to BigInt
354224848179261915075
```
"""
function fibonacciF(n::T; msg=true) where {T<:Integer}

    o = CamiXon.fibonacciF_array(n; msg)[end]

    return o

end
# ..............................................................................
function fibonacciF_array(nmax::T; msg=true) where {T<:Integer}

    n = Int(nmax)
    nc = 92

    if n ≤ nc
        o = T == Int ? glFn_Int[1:1+n] : glFn_BigInt[1:1+n]
    else
        o = _fn_BigInt(n, nc)
        msg && T == Int && println("Warning: fibonacciF autoconverted to BigInt")
    end

    return o

end
#function fibonacciF1(n::T; msg=true) where {T<:Integer}    # kanweg?
#
#    n = Int(n)
#    nc = 92
#
#    if T == Int
#        o = n > nc ? _fn_BigInt(n, nc)[end] : glFn_Int[1+n]
#        msg && n > nc && println("Warning: fibonacciF autoconverted to BigInt")
#    else
#        o = n > nc ? _fn_BigInt(n, nc)[end] : glFn_BigInt[1+n]
#    end
#
#    return o
#
#end
#function fibonacciF_array1(nmax::T; msg=true) where {T<:Integer}     # kanweg?
#
#    n = Int(nmax)
#    nc = 92#
#
#    if T == Int
#        o = n > nc ? _fn_BigInt(n, nc) : glFn_Int[1:1+n]
#        msg && n > nc && println("Warning: fibonacciF autoconverted to BigInt")
#    else
#        o = n > nc ? _fn_BigInt(n, nc) : glFn_BigInt[1:1+n]
#    end
#
#    return o
#
#end

# =================================== harmonic number(n;T) ===============

global glHn_Int = Vector{Rational{Int64}}[
    [1 // 1, 3 // 2, 11 // 6, 25 // 12, 137 // 60, 49 // 20, 363 // 140, 761 // 280, 7129 // 2520, 7381 // 2520,
        83711 // 27720, 86021 // 27720, 1145993 // 360360, 1171733 // 360360, 1195757 // 360360,
        2436559 // 720720, 42142223 // 12252240, 14274301 // 4084080, 275295799 // 77597520,
        55835135 // 15519504, 18858053 // 5173168, 19093197 // 5173168, 444316699 // 118982864,
        1347822955 // 356948592, 34052522467 // 8923714800, 34395742267 // 8923714800,
        312536252003 // 80313433200, 315404588903 // 80313433200, 9227046511387 // 2329089562800,
        9304682830147 // 2329089562800, 290774257297357 // 72201776446800,
        586061125622639 // 144403552893600, 53676090078349 // 13127595717600,
        54062195834749 // 13127595717600, 54437269998109 // 13127595717600,
        54801925434709 // 13127595717600, 2040798836801833 // 485721041551200,
        2053580969474233 // 485721041551200, 2066035355155033 // 485721041551200,
        2078178381193813 // 485721041551200, 85691034670497533 // 19914562703599200,
        12309312989335019 // 2844937529085600, 532145396070491417 // 122332313750680800,
        5884182435213075787 // 1345655451257488800, 5914085889685464427 // 1345655451257488800,
        5943339269060627227 // 1345655451257488800],
    [1 // 1, 5 // 4, 49 // 36, 205 // 144, 5269 // 3600, 5369 // 3600, 266681 // 176400, 1077749 // 705600,
        9778141 // 6350400, 1968329 // 1270080, 239437889 // 153679680, 240505109 // 153679680,
        40799043101 // 25971865920, 40931552621 // 25971865920, 205234915681 // 129859329600,
        822968714749 // 519437318400, 238357395880861 // 150117385017600, 238820721143261 // 150117385017600,
        86364397717734821 // 54192375991353600, 17299975731542641 // 10838475198270720,
        353562301485889 // 221193371393280, 354019312583809 // 221193371393280, 187497409728228241 // 117011293467045120,
        187700554334941861 // 117011293467045120],
    [1 // 1, 9 // 8, 251 // 216, 2035 // 1728, 256103 // 216000, 28567 // 24000, 9822481 // 8232000,
        78708473 // 65856000, 19148110939 // 16003008000, 19164113947 // 16003008000, 25523438671457 // 21300003648000,
        25535765062457 // 21300003648000, 56123375845866029 // 46796108014656000,
        56140429821090029 // 46796108014656000, 56154295334575853 // 46796108014656000,
        449325761325072949 // 374368864117248000],
    [1 // 1, 17 // 16, 1393 // 1296, 22369 // 20736, 14001361 // 12960000, 14011361 // 12960000, 33654237761 // 31116960000,
        538589354801 // 497871360000, 43631884298881 // 40327580160000, 43635917056897 // 40327580160000,
        638913789210188977 // 590436101122560000, 638942263173398977 // 590436101122560000],
    [1 // 1, 33 // 32, 8051 // 7776, 257875 // 248832, 806108207 // 777600000, 268736069 // 259200000,
        4516906311683 // 4356374400000, 144545256245731 // 139403980800000, 105375212839937899 // 101625502003200000,
        105376229094957931 // 101625502003200000],
    [1 // 1, 65 // 64, 47449 // 46656, 3037465 // 2985984, 47463376609 // 46656000000, 47464376609 // 46656000000,
        5584183099672241 // 5489031744000000, 357389058474664049 // 351298031616000000],
    [1 // 1, 129 // 128, 282251 // 279936, 36130315 // 35831808, 2822716691183 // 2799360000000, 940908897061 // 933120000000,
        774879868932307123 // 768464444160000000],
    [1 // 1, 257 // 256, 1686433 // 1679616, 431733409 // 429981696, 168646292872321 // 167961600000000,
        168646392872321 // 167961600000000],
    [1 // 1, 513 // 512, 10097891 // 10077696, 5170139875 // 5159780352, 10097934603139727 // 10077696000000000,
        373997614931101 // 373248000000000],
    [1 // 1, 1025 // 1024, 60526249 // 60466176, 61978938025 // 61917364224, 605263128567754849 // 604661760000000000,
        605263138567754849 // 604661760000000000]
]

# ..............................................................................
global glHn_BigInt = convert(Vector{Vector{Rational{BigInt}}}, glHn_Int)

# ..............................................................................
function _hn_BigInt(n::Int, nc::Int)

    one = big(1)

    o = glHn_BigInt[1][1:nc]
    for m = nc+1:n
        a = o[m-1] + one // big(m)
        Base.push!(o, a)
    end

    return o

end

# ..............................................................................
@doc raw"""
    harmonicNumber(n::T [; msg=true]) where {T<:Integer} 
    harmonicNumber_array(nmax::T [; msg=true]) where {T<:Integer} 

Sum of the reciprocals of the first ``n`` natural numbers
```math
    H_n=\sum_{k=1}^{n}\frac{1}{k}.
```
Integer-overflow protection: for `n > 46` the output is autoconverted to Rational{BigInt}.
By default the capture message is activated: 
"Warning: harmonicNumber autoconverted to Rational{BigInt}". 
### Examples:
```
julia> o = harmonicNumber_array(9); println(o)
Rational{Int64}[1//1, 3//2, 11//6, 25//12, 137//60, 49//20, 363//140, 761//280, 7129//2520]

julia> o = [harmonicNumber(46; msg=true)]; println(o)
Rational{Int64}[5943339269060627227//1345655451257488800]

julia> o = [harmonicNumber(47; msg=true)]; println(o)
Warning: harmonicNumber autoconverted to Rational{BigInt}
Rational{BigInt}[282057509927739620069//63245806209101973600]

julia> harmonicNumber(12) == harmonicNumber(12, 1)
true
```
"""
function harmonicNumber(n::T; msg=true) where {T<:Integer}

    o = CamiXon.harmonicNumber_array(n; msg)[end]

    return o

end

# ..............................................................................
function harmonicNumber_array(nmax::T; msg=true) where {T<:Integer}

    n = Int(nmax)
    nc = 46

    if n ≤ nc
        o = T == Int ? glHn_Int[1][1:n] : glHn_BigInt[1][1:n]
    else
        o = _hn_BigInt(n, nc)
        msg && T == Int && println("Warning: harmonicNumber autoconverted to Rational{BigInt}")
    end

    return o

end
#function harmonicNumber1(n::T; msg=true) where {T<:Integer}
#
#    n = Int(n)
#    nc = 46
#
#    if T == Int
#        o = n > nc ? _hn_BigInt(n, nc)[end] : glHn_Int[1][n]
#        msg && n > nc && println("Warning: harmonicNumber(n) autoconverted to Rational{BigInt}")
#    else
#        o = n > nc ? _hn_BigInt(n, nc)[end] : glHn_BigInt[1][n]
#    end
#
#
#    return o
#
#end
# ..............................................................................
#function harmonicNumber_array1(nmax::T; msg=true) where {T<:Integer}
#
#    n = Int(nmax)
#    nc = 46
#
#    if T == Int
#        o = n > nc ? _hn_BigInt(n, nc) : glHn_Int[1][1:n]
#        msg && n > nc && println("Warning: harmonicNumber(n) autoconverted to Rational{BigInt}")
#    else
#        o = n > nc ? _hn_BigInt(n, nc) : glHn_BigInt[1][1:n]
#    end
#
#    return o
#
#end

# ======================= harmonic number(n, p [; msg=false]) ==========================

function Hn_Int(p::Int, nc::Int)

    o = Rational{Int}[]
    if p > 10
        b = 0 // 1
        for n = 1:nc
            a = 1
            for i = 1:p
                a *= n
            end
            b += 1 // a
            Base.push!(o, b)
        end
    else
        o = glHn_Int[p]
    end

    return o

end
function Hn_BigInt(p::Int, nc::Int)

    nul = big(0)
    one = big(1)

    o = Rational{BigInt}[]
    if p > 10
        b = nul // one
        for k = 1:n
            a = one
            for i = 1:p
                a *= big(k)
            end
            b += one // a
            Base.push!(o, b)
        end
    else
        o = glHn_BigInt[p]
    end

    return o

end
# .......................................................................................
function _hn_BigInt(n::Int, nc::Int, p::Int)

    nul = big(0)
    one = big(1)

    o = CamiXon.Hn_BigInt(p, nc)[1:nc]

    b = nul // one
    for m = 1:n
        a = one
        for i = 1:p
            a *= big(m)
        end
        b += one // a
        Base.push!(o, b)
    end

    return o

end
@doc raw"""
    harmonicNumber(n::T, p::Int [; msg=true]) where {T<:Integer}

Sum of the ``p_{th}`` power of reciprocals of the first ``n`` numbers
```math
    H_{n,p}=\sum_{k=1}^{n}\frac{1}{k^p}.
```
Integer-overflow protection: the output is autoconverted to Rational{BigInt} when required.
By default the capture message is activated: 
"Warning: harmonicNumber autoconverted to Rational{BigInt}". 
### Examples:
```
julia> o = [harmonicNumber(46,1; msg=true)]; println(o)
Rational{Int64}[5943339269060627227//1345655451257488800]

julia> o = [harmonicNumber(47,1; msg=true)]; println(o)
Warning: harmonicNumber autoconverted to Rational{BigInt}"
Rational{BigInt}[280682601097106968469//63245806209101973600]

julia> o = [harmonicNumber6(47,1)]; println(o)
Rational{BigInt}[280682601097106968469//63245806209101973600]

harmonicNumber(12, -3) == faulhaber_summation(12, 3)
  true
```
"""
function harmonicNumber(n::T, p::Int; msg=true) where {T<:Integer}

    n ≠ 0 || return T(0)
    p ≠ 0 || return n

    n = Int(n)
    nc = p < 11 ? length(glHn_Int[p]) : p < 18 ? 4 : p < 25 ? 3 : 0

    if p > 0
        if n ≤ nc
            o = T == Int ? glHn_Int[p][n] : glHn_BigInt[p][n]
        else
            o = _hn_BigInt(n, nc, p)[end]
            msg && T == Int && println("Warning: harmonicNumber autoconverted to Rational{BigInt}")
        end
 #       if T == Int
 #           o = n > nc ? _hn_BigInt(n, nc, p)[end] : Hn_Int(p, nc)[n]
 #           msg && n > nc && println("Warning: harmonicNumber(n, p) autoconverted to Rational{BigInt}")
 #       else
 #           o = n > nc ? _hn_BigInt(n, nc, p)[end] : Hn_BigInt(p, nc)[n]
 #       end
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
julia> triangle_coefficient(3, 4, 5)
1//180180

julia> triangle_coefficient(1//2, 1, 1.5)
1//12
```
"""
function triangle_coefficient(a::Real, b::Real, c::Real)

    (a,b,c) = promote(a,b,c)

    isinteger(a + b + c) || return 0

    A = Int(a + b - c)
    B = Int(b + c - a)
    C = Int(c + a - b)

    A = A ≥ 0 ? bigfactorial(A) : return 0
    B = B ≥ 0 ? bigfactorial(B) : return 0
    C = C ≥ 0 ? bigfactorial(C) : return 0

    num = A * B * C
    den = bigfactorial(Int(a+b+c+1))

    return num//den

end

# ============================ istriangle(a, b, c) =============================

@doc raw"""
    istriangle(a::Real, b::Real, c::Real)

Triangle condition for a triangle of sides `a`, `b` and `c`.

#### Example:
```
julia> istriangle(3, 4, 5)
true

julia> istriangle(1//2, 1, 1.5)
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
julia> p = 5;
julia> texp(1.0, 0.0, 5)
2.7166666666666663

julia> texp(1, 0, 5)
163//60
```
"""
function texp(x::T, a::T, p::Int) where T <: Real

    x = x - a

    V = typeof(x)

    return  V <: Rational ? _texp_int(x, p) : V <: Integer ? _texp_int(x, p) : _texp_real(x, p)

end
