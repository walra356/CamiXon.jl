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
\end{array}\right]=\sum_{j=0}^{k}c_{k-j}^{k}f[n-k+j],
```
for use with the elements of an analytic function, ``f``, tabulated in *forward order*, ``f[n-k], ...,f[n]``  (coefficients in backward order).
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
``k^{th}``*-order lagrangian extrapolation* of the tabulated analytic function ``f(n+x)``
at offset position ``x`` with respect to the position ``n``. Interpolation corresponds to the interval ``-k\le\ x\le 0``.
```math
f[n+x] =\sum_{p=0}^{k}l_p(x)\nabla^pf[n],
```
where ``l_0\equiv 1`` and ``l_p(x) = x(x+1)(x+2)\cdots(x+p-1)/p!``.
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
    summation_range(n, i, k, μ)

Summation range for interpolation position ``i/(μ+1)`` used in ``k^{th}``*-order
lagrangian interpolation* of the anaytic function ``f`` tabulated in forward
order on a uniform grid of ``n`` points, ``f[1],\ \ldots,\ f[n]``;
μ is the interpolation constant.

#### Examples:
```
n = 7; k = 2; μ = 0; m = μ + 1
o = [summation_range(n,i,k,μ) for i=0:(n-1)*m]; println(o)
 UnitRange{Int64}[1:3, 2:4, 3:5, 4:6, 5:7, 5:7, 5:7]
```
"""
function summation_range(n::Int, i::Int, k::Int, μ::Int)
# ================================================================================================
#   summation range for point position i lagrangian interpolation
# ================================================================================================
      m = μ + 1
      0 ≤ i ≤ n*m || error("Error: position index i outside index range 0 ≤ i ≤ n⋅m")
     return i < (n-1-k)*m  ? UnitRange(i÷m+1,i÷m+k+1) : UnitRange(n-k,n)
end

# ==============================================================================

@doc raw"""
    f_diff_function_sequences(f, k::Int, μ=0)

Finite-difference summation sequences of function values given in forward order
for use in ``k^{th}``*-order lagrangian interpolation* of the anaytic function
``f`` tabulated in forward order on a regular grid of ``n`` points, ``f[1], ...,f[n]``;
μ is the interpolation constant. Each sequence consists of ``k⋅m+1`` function values, with ``m=μ+1``.
#### Example:
```
f = [0,1,2,3,4,5,6]
k = 2
o = f_diff_function_sequences(f, k); println(o)
 [[0, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6], [4, 5, 6], [4, 5, 6]]
```
"""
function f_diff_function_sequences(f, k::Int, μ=0)
# ================================================================================================
#   finite-difference function values for interpolation range of lagrangian interpolation
# ================================================================================================
    n = length(f)
    m = μ +1
    return [f[summation_range(n,i,k,μ)] for i=0:(n-1)*m]
end

# ==============================================================================

@doc raw"""
    lagrangian_interpolation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=1, μ=0)

``k^{th}``*-order lagrangian interpolation* of the analytic function ``f`` tabulated
in forward order on a regular grid of ``n`` points, ``f[1],\ \ldots,\ f[n]``;
μ is the interpolation constant (number of intermediate points on the grid).
#### Example:
```
f = [0.0,1,2,3,4,5,6,7]
domain = 0.0..1.0
(X,Y) = lagrangian_interpolation(f, domain; k=2, μ=1); println((X,Y))
 (0.0:0.07142857142857142:1.0, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0])
```
"""
function lagrange_interpolation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=3, μ=0)
# ======================================================================================
#   lagrangian (k+1)-point interpolation at i interpolation points
# ======================================================================================
    n = length(f)
    m = μ + 1

    ∇ = f_diff_weights_array(k)
    l = [f_diff_expansion_coeffs_lagrange(k, x) for x=-k:1/m:0]
    w = [f_diff_expansion_weights(l[i], ∇) for i ∈ eachindex(l)]
    w1 = append!(repeat(w[1:m],n-k-1),w)
    w2 = f_diff_function_sequences(f, k, μ)

    X = range(domain.left, domain.right, length=(n-1)*m+1)
    Y = [w1[i] ⋅ w2[i] for i ∈ eachindex(w1)]

    return X, Y

end

@doc raw"""
    lagrangian_extrapolation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=1, e=1, μ=0)

 ``k^{th}``*-order lagrangian extrapolation* up to position `n+e` with ``μ``
intermediate points of the analytic function ``f`` tabulated in forward order
at ``n`` points, ``f[1],\ \ldots,\ f[n]``.
#### Example:
```
f = [0.0,1,2,3,4,5,6,7]
domain = 0.0..1.0
(X,Y) = lagrangian_extrapolation(f, domain; k=2, e=1, μ=1); println((X,Y))
 (0.0:0.07142857142857142:1.0, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0])
```
"""
function lagrange_extrapolation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=1, e=1, μ=0)
# ======================================================================================
#   lagrangian (k+1)-point interpolation at μ interpolation points
# ======================================================================================
    n = length(f)
    m = μ + 1

    ∇ = f_diff_weights_array(k)
    l = [f_diff_expansion_coeffs_lagrange(k, x) for x=0:1/m:e]
    w1 = [f_diff_expansion_weights(l[i], ∇) for i ∈ eachindex(l)]
    w2 = f_diff_function_sequences(f, k, μ)[end]

    ΔX = (domain.right - domain.left)/((n-1)*m)
    X = range(domain.right, domain.right + ΔX * m*e, length=m*e+1)
    Y = [w1[i] ⋅ w2 for i ∈ eachindex(w1)]

    return X, Y

end

# ===================================== f_diff_expansion_coeffs_differentiation(k, x) ====

@doc raw"""
    f_diff_expansion_coeffs_differentiation(k::Int, x::T) where T<:Real

Finite-difference expansion coefficient vector ``[l_0(x),\ \ldots,\ l_p(x)]`` defining
``k^{th}``*-order lagrangian differentiation* of the tabulated analytic function ``f(n+x)``
at position ``x``.
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
    a = append!([0.0], [1.0/i for i=1:k])
    b = f_diff_expansion_coeffs_lagrange(k, x)

    return polynom_multiplication_coeffs(a, b)[1:k+1]

end

# ================================f_diff_expansion_coeffs_array_differentiation(k, m) ====

@doc raw"""
    lagrange_differentiation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=1, i=0)

``k^{th}``*-order lagrangian differentiation* with ``i`` intermediate points of the function ``f``
tabulated in forward order at ``n`` points, ``f[1],\ \ldots,\ f[n]``.
#### Example:
```
f = [0.0,1,2,3,4,5]
domain = 0.0..5.0
X,Y = lagrangian_differentiation(f, domain; k=2, i = 0); println(X,Y)
 (0.0:1.0:5.0, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
```
"""
function lagrange_differentiation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=1, μ=0)
# ======================================================================================
#   lagrangian (k+1)-point differentiation at i interpolation points
# ======================================================================================
    n = length(f)
    m = μ + 1

    ∇ = f_diff_weights_array(k)
    l = [f_diff_expansion_coeffs_differentiation(k, x) for x=-k:1/m:0]
    w = [f_diff_expansion_weights(l[i], ∇) for i ∈ eachindex(l)]
    w1 = append!(repeat(w[1:m],n-k-1),w)
    w2 = f_diff_function_sequences(f, k, μ)

    X = range(domain.left, domain.right, length=(n-1)*m+1)
    Y = (n-1)/(domain.right-domain.left) .*  [w1[i] ⋅ w2[i] for i ∈ eachindex(w1)]

    return X, Y

end

# ========================== f_diff_expansion_coeffs_adams_moulton(k) ===========


@doc raw"""
    f_diff_expansion_coeffs_adams_moulton(k::Int)

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

    k > 17 ? throw("Error: integer overflow in calculating the 18th-order expansion coefficients (k > 17)") : false

    o = ones(Rational,k+1)
    c = [1//(i+1) for i=1:k]
    for n=1:k
        p = CamiXon.integer_partitions(n)
        a = [1//1 for i ∈ eachindex(p)]
        for i ∈ eachindex(p)
            sgn = sign((-1)^length(p[i]))//1
            a[i] = sgn * permutations_unique_count(p,i)
            for j ∈ eachindex(p[i])
                a[i] *= c[p[i][j]]
            end
        end
        o[n+1] = sum(a)
    end

    return o # Note that D = denominator(gcd(b))

end
