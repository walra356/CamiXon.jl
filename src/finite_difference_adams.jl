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
    fdiff_adams_moulton_expansion_coeff(n::T; msg=true) where {T<:Integer}
    fdiff_adams_moulton_expansion_coeffs(n::T; msg=true) where {T<:Integer}

Finite difference expansion coefficient vector ``β ≡ [β_0(x),\ ⋯,\ β_p(x)]``
defining ``k^{th}``-order Adams-Moulton expansion,

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

julia> o = convert(Vector{Int},(β .* D)); println(o)
[1440, -720, -120, -60, -38, -27]

julia> k=20;
julia> fdiff_adams_moulton_expansion_coeff(k)
Integer-overflow protection: output converted to BigInt
-12365722323469980029//4817145976189747200000
```
"""
function fdiff_adams_moulton_expansion_coeffs(kmax::Int; T=Int, msg=true)

    n = kmax
    nc = 19

    if n ≤ nc
        o = T == Int ? glAMe_Int[1:1+n] : glAMe_BigInt[1:1+n]
    else
        o = _ame_BigInt(n)
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

``k^{th}``-order Adams-Moulton weights vector,
```math
y[n+1] = y[n] + \frac{1}{D}\sum_{j=0}^{k}a^k[j]f[n+1-k+j]
```
The weights are stored in the vector ``a^k \equiv[a_k^k/D,⋯\ a_0^k/D]``
under the convention ``a^k[j] \equiv a_{k-j}^k/D``, where ``a_j^k`` are the
Adams-Moulton weight coefficients and ``D`` the corresponding Adams-Moulton
divisor. By default the output is in Float64, optionally the output is rational,
 with or without specification of the gcd devisor.
#### Example:
```
[create_adams_moulton_weights(k; rationalize=true, devisor=true, T=Int) for k=1:8]
8-element Vector{Tuple{Int64, Int64, Vector{Int64}}}:
 (1, 2, [1, 1])
 (2, 12, [-1, 8, 5])
 (3, 24, [1, -5, 19, 9])
 (4, 720, [-19, 106, -264, 646, 251])
 (5, 1440, [27, -173, 482, -798, 1427, 475])
 (6, 60480, [-863, 6312, -20211, 37504, -46461, 65112, 19087])
 (7, 120960, [1375, -11351, 41499, -88547, 123133, -121797, 139849, 36799])
 (8, 3628800, [-33953, 312874, -1291214, 3146338, -5033120, 5595358, -4604594, 4467094, 1070017])
```
"""
function create_adams_moulton_weights(k::Int; rationalize=false, devisor=false, T=Int)

    β = CamiXon.fdiff_adams_moulton_expansion_coeffs(T(k))

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
    fdiff_adams_bashford_expansion_coeffs(kmax::Int [; T=Int [, msg=true]])
    fdiff_adams_bashford_expansion_coeff(k::Int [; T=Int [, msg=true]])

``(k+1)``-point Adams-Bashford expansion coefficients ``B_p``.

```math
-\frac{∇}{(1-∇)ln(1-∇)}=\sum_{p=0}^{\infty}B_p∇^p=1+\ \frac{1}{2}∇+\ \frac{5}{12}∇^2+\ ⋯.
```
The weights are stored in *forward* order: ``[B_0^k,⋯\ B_k^k]`` -
order of use in summation.
#### Examples:
```
julia> o = fdiff_adams_moulton_expansion_coeffs(5); println(o)
Rational{Int64}[1, -1//2, -1//12, -1//24, -19//720, -3//160]

julia> fdiff_adams_moulton_expansion_coeff(0)
1//1

julia> fdiff_adams_moulton_expansion_coeff(5)
-3//160

julia> fdiff_adams_moulton_expansion_coeff(20)
Integer-overflow protection: output converted to BigInt
-12365722323469980029//4817145976189747200000
```
"""
function fdiff_adams_bashford_expansion_coeffs(nmax::Int; T=Int, msg=true)
# ==============================================================================
#   Adams-Bashford expansion coefficients
# ==============================================================================
    n = Int(nmax)
    nc = 17

    if n ≤ nc
        o = T == Int ? glABe_Int[1:1+n] : glABe_BigInt[1:1+n]
    else
        o = _abe_BigInt(n)
        msg && T == Int && println("Integer-overflow protection: output converted to BigInt")
    end

    return o  # Note that D = denominator(gcd(o))

    end
function fdiff_adams_bashford_expansion_coeff(n::Int; T=Int, msg=true)
    # ==============================================================================
    #   Adams-Bashford expansion coefficients
    # ==============================================================================
      
    o = fdiff_adams_bashford_expansion_coeffs(n; T, msg)
        
    return o[1+n]  

    end
