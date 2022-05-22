# ========================= fdiff_expansion_coeffs_adams_moulton(k) ===========

@doc raw"""
    fdiff_expansion_coeffs_adams_moulton(k [; T=Int])

``k^{th}``-order Adams-Moulton expansion coefficients,

```math
-\frac{∇}{ln(1-∇)} = \sum_{p=0}^{\infty}b_p∇^p= 1 - \frac{1}{2}∇ - \frac{1}{12}∇^2 - \frac{1}{24}∇^3 +⋯.
```
The weights are stored in *forward* order: ``[b_0^k,⋯\ b_k^k]`` -
order of use in summation.
#### Examples:
```
k = 5
b = fdiff_expansion_coeffs_adams_moulton(k::Int); println(b)
 Rational[1//1, -1//2, -1//12, -1//24, -19//720, -3//160]

D = denominator(gcd(b)); println(D)
 1440

o = convert(Vector{Int},(b .* D)); println(o)
 [1440, -720, -120, -60, -38, -27]
```
"""
function fdiff_expansion_coeffs_adams_moulton(k::Int; T=Int)
# ==============================================================================
#   Adams-Moulton expansion coefficients
# ==============================================================================
    o = Base.zeros(Rational{T},k+1)
    s::Rational{T} = 0//1

    o[1] = 1//1

    b = Base.copy(o)
    a = Base.append!([s],[1//(i+1) for i=1:k])

    for p=1:k
        b = CamiXon.polynom_product_expansion(a, b, k)
        Base.isodd(p) ? o = o .- b : o = o .+ b
    end

    return o  # Note that D = devisor(gcd(o))

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

    ∇ = CamiXon.fdiff_weights_array(k)
    #∇ = fdiffs(k, forward=false)
    β = CamiXon.fdiff_expansion_coeffs_adams_moulton(k; T)

    o = CamiXon.fdiff_expansion_weights(β, ∇; notation=bwd)

    if rationalize
        D = Base.denominator(Base.gcd(o))       # Adams-Moulton devisor
        o = devisor ? (k, D, Base.round.(T, o * D)) : o
    else
        o = Base.convert(Vector{Float64},o)
    end

    return o

end

# ======================== fdiff_expansion_coeffs_adams_bashford(k) ===========

@doc raw"""
    fdiff_expansion_coeffs_adams_bashford(k [; T=Int])

``(k+1)``-point Adams-Bashford expansion coefficients ``B_p``.

```math
-\frac{∇}{(1-∇)ln(1-∇)}=\sum_{p=0}^{\infty}B_p∇^p=1+\ \frac{1}{2}∇+\ \frac{5}{12}∇^2+\ ⋯.
```
The weights are stored in *forward* order: ``[B_0^k,⋯\ B_k^k]`` -
order of use in summation.
#### Examples:
```
k = 5
o = fdiff_expansion_coeffs_adams_bashford(k); println(o)
 Rational{Int64}[1//1, 1//2, 5//12, 3//8, 251//720, 95//288]
```
"""
function fdiff_expansion_coeffs_adams_bashford(k::Int; T=Int)
# ==============================================================================
#   Adams-Bashford expansion coefficients
# ==============================================================================
    a = Base.ones(Rational{T},k+1)

    b = CamiXon.fdiff_expansion_coeffs_adams_moulton(k; T)
    o = CamiXon.polynom_product_expansion(a, b, k)

    return o  # Note that D = denominator(gcd(o))

end
