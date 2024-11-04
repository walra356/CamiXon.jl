# SPDX-License-Identifier: MIT

# author: Jook Walraven - 14-2-2023

# ==============================================================================
#                              finite_difference_adams.jl
# ==============================================================================

# ========================= fdiff_adams_moulton_expansion_coeffs(k) ===========

global glAMe_Int = [1 // 1, -1 // 2, -1 // 12, -1 // 24, -19 // 720, -3 // 160, -863 // 60480, -275 // 24192,
    -33953 // 3628800, -8183 // 1036800, -3250433 // 479001600, -4671 // 788480,
    -13695779093 // 2615348736000, -2224234463 // 475517952000,
    -132282840127 // 31384184832000, -2639651053 // 689762304000,
    -111956703448001 // 32011868528640000, -50188465 // 15613165568,
    -2334028946344463 // 786014494949376000, -301124035185049 // 109285437800448000]

# ..............................................................................
global glAMe_BigInt = convert(Vector{Rational{BigInt}}, glAMe_Int)

# ..............................................................................
function _ame_BigInt(k::Int)

    nul = big(0)
    one = big(1)

    o = Base.zeros(Rational{BigInt}, k + 1)
    o[1] = big(1) // big(1)

    a = [one // big(i + 1) for i = 0:k]
    a[1] = nul // one

    b = Base.copy(o)

    for p = 1:k
        b = CamiMath.polynom_product_expansion(a, b, k)
        Base.isodd(p) ? o = o .- b : o = o .+ b
    end

    return o

end
# ..............................................................................

@doc raw"""
    fdiff_adams_moulton_expansion_coeff(k::Int; T=Int, msg=true)
    fdiff_adams_moulton_expansion_coeffs(k::Int; T=Int, msg=true)

Finite difference expansion coefficient vector ``β ≡ [β_0(x),\ ⋯,\ β_k(x)]``.
Note the *forward* vector ordering, which is the order of use in the summation below,

```math
-\frac{∇}{ln(1-∇)}
= \sum_{p=0}^{\infty}β_p∇^p
= 1 - \frac{1}{2}∇ - \frac{1}{12}∇^2 - \frac{1}{24}∇^3 +⋯.
```
#### Examples:
```
julia> k = 5;
julia> β = fdiff_adams_moulton_expansion_coeffs(k); println(β)
Rational{Int64}[1//1, -1//2, -1//12, -1//24, -19//720, -3//160]

julia> fdiff_adams_moulton_expansion_coeff(k)
-3//160

julia> D = denominator(gcd(β))
1440

julia> convert(Vector{Int},(β .* D)); println(o)
[1440, -720, -120, -60, -38, -27]

julia> k = 20;
julia> fdiff_adams_moulton_expansion_coeff(k)
Integer-overflow protection: output converted to BigInt
-12365722323469980029//4817145976189747200000
```
"""
function fdiff_adams_moulton_expansion_coeffs(k::Int; T=Int, msg=true)

    nc = 19

    if k ≤ nc
        o = T == Int ? glAMe_Int[1:1+k] : glAMe_BigInt[1:1+k]
    else
        o = _ame_BigInt(k)
        msg && T == Int && println("Integer-overflow protection: output converted to BigInt")
    end

    return o

end
function fdiff_adams_moulton_expansion_coeff(k::Int; T=Int, msg=true)

    o = fdiff_adams_moulton_expansion_coeffs(k; T, msg)

    return o[1+k]

end

# ========================== create_adams_moulton_weights(k)====================

@doc raw"""
    create_adams_moulton_weights(k::Int [; rationalize=false [, devisor=false [, T=Int]]])

``k^{th}``-order Adams-Moulton weights vector ``a^k \equiv[a_k^k,⋯\ a_0^k]``.  
Note the *reversed* vector ordering, which is the order of use in the summation below,

```math
y[n+1] = y[n] + \frac{1}{D}\sum_{j=0}^{k}a^k[j]f[n+1-k+j],
```
where ``a^k[j] \equiv a_{k-j}^k``. The ``a_j^k`` are the
Adams-Moulton weight coefficients and ``D`` is the corresponding Adams-Moulton
divisor. By default the output is in Float64, optionally the output is rational,
 with or without specification of the gcd devisor.
#### Example:
```
julia> [create_adams_moulton_weights(k; rationalize=true, devisor=true, T=Int) for k=1:5]
8-element Vector{Tuple{Int64, Int64, Vector{Int64}}}:
 (1, 2, [1, 1])
 (2, 12, [-1, 8, 5])
 (3, 24, [1, -5, 19, 9])
 (4, 720, [-19, 106, -264, 646, 251])
 (5, 1440, [27, -173, 482, -798, 1427, 475])

julia> k = 5;
julia> w = create_adams_moulton_weights(k; rationalize=true, devisor=true); println(w)
(5, 1440, [27, -173, 482, -798, 1427, 475])

julia> w = create_adams_moulton_weights(k; rationalize=true, devisor=false); println(w)
Rational{Int64}[3//160, -173//1440, 241//720, -133//240, 1427//1440, 95//288]

julia> w = create_adams_moulton_weights(k; rationalize=false); println(w)
[0.01875, -0.12013888888888889, 0.3347222222222222, -0.5541666666666667, 0.9909722222222223, 0.3298611111111111]
```
"""
function create_adams_moulton_weights(k::Int; rationalize=false, devisor=false, T=Int)

    β = CamiXon.fdiff_adams_moulton_expansion_coeffs(k; T)

    o = fdiff_expansion_weights(β)

    if rationalize
        D = Base.denominator(Base.gcd(o))       # Adams-Moulton devisor
        o = devisor ? (k, D, Base.round.(T, o * D)) : o
    else
        o = Base.convert(Vector{Float64},o)
    end

    return o

end

# ======================== fdiff_adams_bashford_expansion_coeffs(k) ===========

global glABe_Int = Rational{Int64}[1, 1//2, 5//12, 3//8, 251//720, 95//288, 19087//60480, 
    5257//17280, 1070017//3628800, 25713//89600, 26842253//95800320, 4777223//17418240, 
    703604254357//2615348736000, 106364763817//402361344000, 1166309819657//4483454976000, 
    25221445//98402304, 8092989203533249//32011868528640000, 85455477715379//342372925440000]

# ..............................................................................
global glABe_BigInt = convert(Vector{Rational{BigInt}}, glABe_Int)

# ..............................................................................
function _abe_BigInt(k)

    a = Base.ones(Rational{BigInt},1+k)
    
    b = CamiXon.fdiff_adams_moulton_expansion_coeffs(k; msg=false)
    o = CamiMath.polynom_product_expansion(a, b, k)

    return o  # Note that D = denominator(gcd(o))

end
# ..............................................................................

@doc raw"""
    fdiff_adams_bashford_expansion_coeff(k::Int [; T=Int [, msg=true]])
    fdiff_adams_bashford_expansion_coeffs(k::Int [; T=Int [, msg=true]])

``(k+1)``-point Adams-Bashford expansion coefficients ``B_k \equiv [B_0^k,⋯\ B_k^k]``. 
Note the *forward* vector ordering, which is the order of use in the summation below,

```math
-\frac{∇}{(1-∇)ln(1-∇)}=\sum_{p=0}^{\infty}B_p∇^p=1+\ \frac{1}{2}∇+\ \frac{5}{12}∇^2+\ ⋯.
```
#### Examples:
```
julia> o = fdiff_adams_bashford_expansion_coeffs(5); println(o)
Rational{Int64}[1, 1//2, 5//12, 3//8, 251//720, 95//288]

julia> fdiff_adams_bashford_expansion_coeff(0)
1//1

julia> fdiff_adams_bashford_expansion_coeff(5)
95//288

julia> fdiff_adams_bashford_expansion_coeff(20)
Integer-overflow protection: output converted to BigInt
8136836498467582599787//33720021833328230400000
```
"""
function fdiff_adams_bashford_expansion_coeffs(k::Int; T=Int, msg=true)
# ==============================================================================
#   Adams-Bashford expansion coefficients
# ==============================================================================
    nc = 17

    if k ≤ nc
        o = T == Int ? glABe_Int[1:1+k] : glABe_BigInt[1:1+k]
    else
        o = _abe_BigInt(k)
        msg && T == Int && println("Integer-overflow protection: output converted to BigInt")
    end

    return o  # Note that D = denominator(gcd(o))

    end
function fdiff_adams_bashford_expansion_coeff(k::Int; T=Int, msg=true)
      
    o = fdiff_adams_bashford_expansion_coeffs(k; T, msg)
        
    return o[1+k]  

    end

    # ========================== create_adams_bashford_weights(k)====================

@doc raw"""
    create_adams_bashford_weights(k::Int [; rationalize=false [, devisor=false [, T=Int]]])

``k^{th}``-order Adams-Bashford weights vector ``b^k \equiv[b_k^k,⋯\ b_0^k]``. 
Note the *reversed* order, which corresponds to the order of use in the summation below,

```math
y[n+1] = y[n] + \frac{1}{D}\sum_{j=0}^{k}b^k[j]f[n+1-k+j],
```
where ``b^k[j] \equiv b_{k-j}^k``. The ``b_j^k`` are the
Adams-Bashford weight coefficients, with ``D`` the corresponding Adams-Moulton
divisor. By default the output is in Float64, optionally the output is rational,
with or without specification of the gcd devisor.
#### Example:
```
julia> [create_adams_bashford_weights(k; rationalize=true, devisor=true, T=Int) for k=1:5]
8-element Vector{Tuple{Int64, Int64, Vector{Int64}}}:
 (1, 2, [-1, 3])
 (2, 12, [5, -16, 23])
 (3, 24, [-9, 37, -59, 55])
 (4, 720, [251, -1274, 2616, -2774, 1901])
 (5, 1440, [-475, 2877, -7298, 9982, -7923, 4277])

julia> k = 5;

julia> w = create_adams_bashford_weights(k; rationalize=true, devisor=true); println(w)
(5, 1440, [-475, 2877, -7298, 9982, -7923, 4277])

julia> w = create_adams_bashford_weights(k; rationalize=true, devisor=false); println(w)
Rational{Int64}[-95//288, 959//480, -3649//720, 4991//720, -2641//480, 4277//1440]

julia> w = create_adams_bashford_weights(k; rationalize=false); println(w)
[-0.3298611111111111, 1.9979166666666666, -5.0680555555555555, 6.9319444444444445, -5.502083333333333, 2.970138888888889]
```
"""
function create_adams_bashford_weights(k::Int; rationalize=false, devisor=false, T=Int)

B = CamiXon.fdiff_adams_bashford_expansion_coeffs(k; T)

o = fdiff_expansion_weights(B)

if rationalize
    D = Base.denominator(Base.gcd(o))       # Adams-Bashford devisor = Adams-Moulton devisor
    o = devisor ? (k, D, Base.round.(T, o * D)) : o
else
    o = Base.convert(Vector{Float64},o)
end

return o

end