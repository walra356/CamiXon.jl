@doc raw"""
    f_diff_weight(k, j)

Weight coefficient
```math
c_{j}^{k}=(-1)^{j}\binom{k}{j},
```
of the ``k^{th}``-order finite difference operator ``\nabla^k`` and corresponding to the function value ``f[n-j]``.
#### Example:
```
k = 5; j = 3
f_diff_weight(k, j)
 -10
```
"""
f_diff_weight(k::Int, j::Int) = iseven(j) ? Base.binomial(k,j) : -Base.binomial(k,j)

# ==============================================================================

@doc raw"""
    f_diff_weights(k)

Weight vector ``[c_k^k,\ \ldots,\ c_0^k]`` defining the ``k^{th}``-order finite difference operator,
```math
\nabla^{k}f[n]	=[c_{k}^{k},\thinspace c_{k-1}^{k},\thinspace\ldots,c_{0}^{k}]\left[\begin{array}{c}
f[n-k]\\
\vdots\\
f[n]
\end{array}\right]=\sum_{j=0}^{k}c_{k-j}^{k}f[n-k+j].
```
This form is recommended for use with any analytic function, ``f``, tabulated in *forward order*, ``f[n-k], ...,f[n]``  (coefficients in backward order).
#### Example:
```
k = 3
f_diff_weights(k)
4-element Vector{Int64}:
  1
 -3
  3
 -1
```
"""
f_diff_weights(k::Int) = [f_diff_weight(k, k-j) for j=0:k]

# ==============================================================================

@doc raw"""
    f_diff_weights_array(kmax)

Collection of weight vectors ``[c_k^k,\ \ldots,\ c_0^k]``  defining the finite difference operators ``\nabla^0,\ \ldots,\ \nabla^k``.
#### Example:
```
kmax = 3
∇ = f_diff_weights_array(kmax)
4-element Vector{Vector{Int64}}:
 [1]
 [-1, 1]
 [1, -2, 1]
 [-1, 3, -3, 1]
```
"""
f_diff_weights_array(kmax::Int) = [f_diff_weights(k)  for k=0:kmax]

# ==============================================================================

@doc raw"""
    f_diff_expansion_weights(a, ∇)

Weight vector ``[b_k^k,\ ,\ldots,\ b_0^k]`` corresponding to the expansion coefficients
``[a_0^k,\ ,\ldots,\ a_k^k]`` of the ``k^{th}``-order finite-difference expansion,

```math
\sum_{p=0}^{k}a_{p}\nabla^{p}f[n]=\sum_{j=0}^{k}b_{k-j}^{k}f[n-k+j],
```

where ``f[n-k], ...,f[n]`` are elements of the analytic function ``f`` tabulated in *forward order*.
Note the difference in ordering between the finite-difference expansion *coefficients*,
``a_{0},\ \ldots,\ a_{k}``, and the finite-difference expansion *weights*, ``b_k^{k},\ \ldots,\ b_0^{k}``.
Note further the difference in ``k`` dependence: the *weights*,
``b_j^k``, are ``k``*-dependent*, whereas the *coefficients*, ``a_j``, are not.
#### Example:
```
k=5
∇ = f_diff_weights_array(k)
a = UnitRange(0,k)
b = f_diff_expansion_weights(a, ∇)
6-element Vector{Int64}:
  -5
  29
 -69
  85
 -55
  15
```
"""
function f_diff_expansion_weights(coeffs, ∇)
# ======================================================================================
#   function weights of finite-difference summation
# ======================================================================================
    k = length(coeffs)-1
    #return [sum([coeffs[1+p] * ∇[1+p][1+p-j] for p=j:k]) for j=0:k]
    return [sum([coeffs[1+p] * ∇[1+p][1+p-j] for p=j:k]) for j=k:-1:0]
end

# ==============================================================================

@doc raw"""
    f_diff_expansion_coeffs_lagrange(k::Int, x::T) where T<:Real

Finite-difference expansion coefficient vector ``[l_0(x),\ \ldots,\ l_p(x)]`` defining
``k^{th}``-order lagrangian interpolation of the tabulated analytic function ``f(n+x)``
at offset position ``x`` with respect to position ``n``,
```math
f[n+x] = (1 - \nabla)^{-x} f[n] \equiv \sum_{p=0}^{\infty}l_p(x)∇^p f[n],
```
where ``l_0\equiv 1`` and ``l_p(x) = x(x+1)(x+2)\cdots(x+p-1)/p!``. Interpolation corresponds to the interval
``-k\le\ x\le 0``; extrapolation to ``x\ge 0``.
#### Examples:
```
k = 5; x = 1
l = f_diff_expansion_coeffs_lagrange(k,x); println(l)
 [1, 1, 1, 1, 1, 1]
```
"""
function f_diff_expansion_coeffs_lagrange(k::Int, x::T) where T<:Real
# ======================================================================================
#   f_difference expansion coefficients for lagrange interpolation for position n+x
# ======================================================================================
    l = ones(T,k+1)
    x == 1 ? l : x ==0 ? (for i=2:k+1; l[i] = 0 end) : (for i=1:k; l[i+1] = l[i]*(x+i-1)/i end)
    return l
end

# ==============================================================================

@doc raw"""
    summation_range(n, i, k, m)

Summation range for interpolation position ``0\le i/m \le 1`` used in ``k^{th}``-order
lagrangian interpolation of the anaytic function ``f`` tabulated in forward
order on a uniform grid of ``n`` points, ``f[1],\ \ldots,\ f[n]``;
``m`` is the multiplier defining the interpolation grid size.

#### Examples:
```
n = 7; k = 2; m = 1
o = [summation_range(n,i,k,m) for i=0:(n-1)*m]; println(o)
 UnitRange{Int64}[1:3, 2:4, 3:5, 4:6, 5:7, 5:7, 5:7]
```
"""
function summation_range(n::Int, i::Int, k::Int, m::Int)
# ================================================================================================
#   summation range for point position i lagrangian interpolation
# ================================================================================================
      0 ≤ i ≤ n*m || error("Error: position index i outside index range 0 ≤ i ≤ n⋅m")
     return i < (n-1-k)*m  ? UnitRange(i÷m+1,i÷m+k+1) : UnitRange(n-k,n)
end

# ==============================================================================

@doc raw"""
    f_diff_function_sequences(f, k::Int, m=1)

Finite-difference summation sequences of function values given in forward order
for use in ``k^{th}``-order lagrangian interpolation of the anaytic function
``f`` tabulated in forward order on a regular grid of ``n`` points, ``f[1], ...,f[n]``;
``m`` is the multiplier defining the interpolation grid size. Each sequence consists of ``k⋅m+1`` function values.
#### Example:
```
f = [0,1,2,3,4,5,6]
k = 2
o = f_diff_function_sequences(f, k); println(o)
 [[0, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6], [4, 5, 6], [4, 5, 6]]
```
"""
function f_diff_function_sequences(f, k::Int, m=1)
# ================================================================================================
#   finite-difference function values for interpolation range of lagrangian interpolation
# ================================================================================================
    n = length(f)
    return [f[summation_range(n,i,k,m)] for i=0:(n-1)*m]
end

# ==============================================================================

@doc raw"""
    lagrangian_interpolation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=1, m=1)

``k^{th}``-order lagrangian *interpolation* of the analytic function ``f`` tabulated
in forward order on a regular grid of ``n`` points, ``f[1],\ \ldots,\ f[n]``;
``m`` is the multiplier defining the interpolation grid size.
#### Example:
```
f = [0.0,1,2,3,4,5,6,7]
domain = 0.0..1.0
(X,Y) = lagrangian_interpolation(f, domain; k=2, m=2); println((X,Y))
 (0.0:0.07142857142857142:1.0, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0])
```
"""
function lagrange_interpolation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=3, m=1)
# ======================================================================================
#   lagrangian (k+1)-point interpolation at i interpolation points
# ======================================================================================
    n = length(f)

    ∇ = f_diff_weights_array(k)
    l = [f_diff_expansion_coeffs_lagrange(k, x) for x=-k:1/m:0]
    w = [f_diff_expansion_weights(l[i], ∇) for i ∈ eachindex(l)]
    w1 = append!(repeat(w[1:m],n-k-1),w)
    w2 = f_diff_function_sequences(f, k, m)

    X = range(domain.left, domain.right, length=(n-1)*m+1)
    Y = [w1[i] ⋅ w2[i] for i ∈ eachindex(w1)]

    return X, Y

end

@doc raw"""
    lagrangian_extrapolation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=1, e=1, m=1)

``k^{th}``-order lagrangian *extrapolation* up to position ``n+e`` of the analytic function
``f`` tabulated in forward order at ``n`` points, ``f[1],\ \ldots,\ f[n]``;
``m`` is the multiplier defining the interpolation grid size.
#### Example:
```
f = [0.0,1,2,3,4,5,6,7]
domain = 0.0..1.0
(X,Y) = lagrangian_extrapolation(f, domain; k=2, e=1, m=2); println((X,Y))
 (0.0:0.07142857142857142:1.0, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0])
```
"""
function lagrange_extrapolation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=1, e=1, m=1)
# ======================================================================================
#   lagrangian (k+1)-point interpolation at μ interpolation points
# ======================================================================================
    n = length(f)

    ∇ = f_diff_weights_array(k)
    l = [f_diff_expansion_coeffs_lagrange(k, x) for x=0:1/m:e]
    w1 = [f_diff_expansion_weights(l[i], ∇) for i ∈ eachindex(l)]
    w2 = f_diff_function_sequences(f, k, m)[end]

    ΔX = (domain.right - domain.left)/((n-1)*m)
    X = range(domain.right, domain.right + ΔX * m*e, length=m*e+1)
    Y = [w1[i] ⋅ w2 for i ∈ eachindex(w1)]

    return X, Y

end

# ===================================== f_diff_expansion_coeffs_differentiation(k, x) ====

@doc raw"""
    f_diff_expansion_coeffs_differentiation(k::Int, x::T) where T<:Real

Finite-difference expansion coefficient vector ``[l_0^{\prime}(x),\ \ldots,\ l_p^{\prime}(x)]`` defining
``k^{th}``-order lagrangian *differentiation* of the tabulated analytic function ``f(n+x)``
at position ``x``,
```math
\frac{df}{dx}[n+x]=\sum_{p=0}^{k}l_{p}^{\prime}(x)\nabla^{p}f[n]
```
#### Example:
```
k = 2; x = 0
o = f_diff_expansion_coeffs_differentiation(k,x); println(o)
 [0.0, 1.0, -1.5]
```
"""
function f_diff_expansion_coeffs_differentiation(k::Int, x::T) where T<:Real
# ======================================================================================
#   finite difference expansion coeffs for differentiation in interval -k ≤ x ≤ 0
# ======================================================================================
    a = prepend!([1//i for i=1:k],[0//1])
    b = f_diff_expansion_coeffs_lagrange(k, x)

    return polynom_multiplication_coeffs(a, b)[1:k+1]

end

# ================================f_diff_expansion_coeffs_array_differentiation(k, m) ====

@doc raw"""
    lagrange_differentiation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=1, m=1)

``k^{th}``-order lagrangian *differentiation* of the analytic function ``f``, tabulated
in forward order on a regular grid of ``n`` points, ``f[1],\ \ldots,\ f[n]``;
``m`` is the multiplier for intermediate positions
#### Example:
```
f = [0.0,1,2,3,4,5]
domain = 0.0..5.0
X,Y = lagrangian_differentiation(f, domain; k=2, i = 0); println(X,Y)
 (0.0:1.0:5.0, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
```
"""
function lagrange_differentiation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=1, m=1)
# ======================================================================================
#   lagrangian (k+1)-point differentiation at i interpolation points
# ======================================================================================
    n = length(f)

    ∇ = f_diff_weights_array(k)
    l = [f_diff_expansion_coeffs_differentiation(k, x) for x=-k:1/m:0]
    w = [f_diff_expansion_weights(l[i], ∇) for i ∈ eachindex(l)]
    w1 = append!(repeat(w[1:m],n-k-1),w)
    w2 = f_diff_function_sequences(f, k, m)

    X = range(domain.left, domain.right, length=(n-1)*m+1)
    Y = (n-1)/(domain.right-domain.left) .*  [w1[i] ⋅ w2[i] for i ∈ eachindex(w1)]

    return X, Y

end

# ========================== f_diff_expansion_coeffs_adams_moulton(k) ===========


@doc raw"""
    f_diff_expansion_coeffs_adams_moulton(k::Int)

Adams-Moulton finite-difference expansion coefficients (restricted to order k < 18),

```math
-\frac{\nabla}{ln(1-\nabla)} = \sum_{p=0}^{\infty}b_p\nabla^p= 1 - \frac{1}{2}\nabla - \frac{1}{12}\nabla^2 - \frac{1}{24}\nabla^3 +\cdots.
```
#### Examples:
```
k = 5
b = f_diff_expansion_coeffs_adams_moulton(k::Int); println(b)
 Rational[1//1, -1//2, -1//12, -1//24, -19//720, -3//160]

D = denominator(gcd(b)); println(D)
 1440

o = convert(Vector{Int},(b .* D)); println(o)
 [1440, -720, -120, -60, -38, -27]
```
"""
function f_diff_expansion_coeffs_adams_moulton(k::Int)
# =====================================================================================
#   Adams-Moulton expansion coefficients limited to order k < 18
# =====================================================================================
    k < 18 || error("Error: integer overflow in calculating the 19th-order expansion coefficients (k > 17)")

    a = Base.prepend!([1//(i+1) for i=1:k],[0//1])
    o = Base.zeros(Rational{Int}, k+1); o[1] = 1//1
    b = Base.copy(o)

    for p=1:k
        b = CamiXon.polynom_multiplication_coeffs(a, b)[1:k+1]
        Base.isodd(p) ? o = o .- b : o = o .+ b
    end

    return o  # Note that D = denominator(gcd(o))

end

# ========================== f_diff_expansion_coeffs_adams_bashford(k) ===========

@doc raw"""
    f_diff_expansion_coeffs_adams_bashford(k::Int)

Adams-Bashford finite-difference-expansion coefficients ``B_p`` (restricted to order k < 18)

```math
-\frac{\nabla}{(1-\nabla)ln(1-\nabla)}=\sum_{p=0}^{\infty}B_p\nabla^p=1+\ \frac{1}{2}∇+\ \frac{5}{12}∇^2+\ \cdots.
```
#### Examples:
```
k = 5
o = f_diff_expansion_coeffs_adams_bashford(k::Int); println(o)
 Rational{Int64}[1//1, 1//2, 5//12, 3//8, 251//720, 95//288]
```
"""
function f_diff_expansion_coeffs_adams_bashford(k::Int)
# =====================================================================================
#   Adams-Bashford expansion coefficients restricted to order k < 18
# =====================================================================================
    k < 18 || error("Error: integer overflow in calculating the 19th-order expansion coefficients (k > 17)")

    a = ones(Int, k+1) #f_diff_expansion_coeffs_lagrange(k,1)
    b = f_diff_expansion_coeffs_adams_moulton(k)
    o = polynom_multiplication_coeffs(a, b)[1:k+1]

    return o  # Note that D = denominator(gcd(o))

end
