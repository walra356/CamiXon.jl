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
where ``f[n-k], ...,f[n]`` are elements of the analytic function ``f`` tabulated in *normal ordering*.
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
\sum_{p=0}^{k}a_{p}\nabla^{p}f[n]=\sum_{j=0}^{k}b_{j}^{k}f[n-j]=\sum_{j=0}^{k}b_{k-j}^{k}f[n-k+j],
```

where ``f[n-k], ...,f[n]`` are elements of the analytic function ``f`` tabulated in *normal ordering*.
Note the difference in ordering between the finite-difference expansion *coefficients*,
``a_{0},\ \ldots,\ a_{k}``, and the finite-difference expansion *weights*, ``b_k^{k},\ \ldots,\ b_0^{k}``.
Note that the *coefficients*, ``a_j``, are *independent of* ``k``, whereas the *weights*,
``b_j^k``, are ``k``*-dependent*.
#### Example:
```
k=5
∇ = f_diff_weights_array(k)
a = UnitRange(0,k)
b = f_diff_expansion_weights(a, ∇)
6-element Vector{Int64}:
  15
 -55
  85
 -69
  29
  -5
```
"""
function f_diff_expansion_weights(coeffs, ∇)
# ======================================================================================
#   function weights of finite-difference summation
# ======================================================================================
    k = length(coeffs)-1
    return [sum([coeffs[1+p] * ∇[1+p][1+p-j] for p=j:k]) for j=0:k]
end

# ==============================================================================

@doc raw"""
    f_diff_expansion_coeffs_interpolation(k::Int, x::T) where T<:Real

Finite-difference expansion coefficient vector ``[l_0(x),\ \ldots\,\ l_p(x)]`` defining
``(k+1)``*-point lagrangian interpolation* of the tabulated analytic function ``f(n+x)``
at offset position ``x`` with respect to the position ``n``, with ``-k\le x\le 0``,
```math
f[n+x] =\sum_{p=0}^{k}l_p(x)\nabla^pf[n],
```
where ``l_0\equiv 1`` and ``l_p(x) = x(x+1)(x+2)\cdots(x+p-1)/p!``.
#### Examples:
```
k=3
∇ = f_diff_weights_array(k)
x=-1
l = f_diff_expansion_coeffs_interpolation(k,x)
r = f_diff_expansion_weights(l, ∇)
println(l,r)
 [1, -2, 1, 0][0, 0, 1, 0]
```
"""
function f_diff_expansion_coeffs_interpolation(k::Int, x::T) where T<:Real
# ======================================================================================
#   f_difference expansion parameters for the interpolation interval -k ≤ x ≤ 0
# ======================================================================================
    x > 0 ? error("Error: outside interpolation range (x > 0)") :
    x < -k ? error("Error: outside interpolation range (x < $(-k))") :
    l = ones(T,k+1)
    for i=1:k
        l[i+1] = l[i]*(-(x+k)+i-1)/i
    end
    return l
end

# ==============================================================================

@doc raw"""
    f_diff_expansion_coeffs_array_interpolation(k::Int, m::Int)

Finite-difference expansion coefficient vector ``[l_0(x),\ \ldots\,\ l_p(x)]`` defining
``(k+1)``*-point lagrangian interpolation* of the tabulated analytic function ``f(n+x)``
for ``k*m+1`` x positions.
####
```
k = 2; m = 2
o = f_diff_expansion_coeffs_array_interpolation(k,m); println(o)
 [[1.0, 0.0, 0.0], [1.0, -0.5, -0.125], [1.0, -1.0, -0.0], [1.0, -1.5, 0.375], [1.0, -2.0, 1.0]]

m = 1
o = f_diff_expansion_coeffs_array_interpolation(k,m); println(o)
 [[1, 0, 0], [1, -1, 0], [1, -2, 1]]
```
"""
function f_diff_expansion_coeffs_array_interpolation(k::Int, m::Int)
    p = m ≠ 1 ? [i/m-k for i=0:k*m] : [i-k for i=0:k]
    return [f_diff_expansion_coeffs_interpolation(k, p[i]) for i ∈ eachindex(p)]
end

# ==============================================================================

@doc raw"""
    summation_ranges(n, k, i, m)

Summation ranges for interpolation positions ``i = 0,\ 1,\ \ldots,\ (n-1)*m`` as used in
``k^{th}``-order lagrangian interpolation of an anaytic function
``f`` tabulated in normal ordering on a uniform grid of ``n`` points, f[1], ...,f[n].
#### Examples:
```@docs
n = 7; k = 2; m = 1
o = [summation_ranges(n,i,k,m) for i=0:(n-1)*m]; println(o)
 UnitRange{Int64}[1:3, 2:4, 3:5, 4:6, 5:7, 5:7, 5:7]
```
"""
function summation_ranges(n::Int, i::Int, k::Int, m::Int)
# ================================================================================================
#   summation range for point position i lagrangian interpolation
# ================================================================================================
    return i < (n-1-k)*m  ? UnitRange(i÷m+1,i÷m+k+1) : UnitRange(n-k,n)
end

# ==============================================================================

@doc raw"""
    f_diff_function_sequences(f, n::Int, k::Int, m::Int)

Finite-difference interpolation sequences (of ``k+1`` function values in normal order) for
``k^{th}``-order lagrangian intepolation of the anaytic function f(i) tabulated in normal order
at n points, f[1], ...,f[n].
#### Example:
```
f = [0,1,2,3,4,5,6]
k = 2; m = 1
o = f_diff_function_sequences(f, k, m); println(o)
 [[0, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6], [4, 5, 6], [4, 5, 6]]
```
"""
function f_diff_function_sequences(f, k::Int, m::Int)
# ================================================================================================
#   finite-difference function values for interpolation range of lagrangian interpolation
# ================================================================================================
    n = length(f)
    return [f[summation_ranges(n,i,k,m)] for i=0:(n-1)*m]
end

# ==============================================================================

@doc raw"""
    f_diff_expansion_weights_array(n::Int, k::Int, m::Int, l::Vector{Vector{T}}) where T<:Real

Weight factors for ``k^{th}-order finite-difference expansion for finite-difference expansion
coefficients ``l`` on uniform grid of ``(n-1)*m+1`` points.
#### Example:
```
n = 7; k = 3; m = 1
l = f_diff_expansion_coeffs_array_interpolation(k, m)
o = f_diff_expansion_weights_array(n, k, m, l); println(o)
 [[1, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
```
"""
function f_diff_expansion_weights_array(n::Int, k::Int, m::Int, l::Vector{Vector{T}}) where T<:Real
# ================================================================================================
#   function weights of finite-difference expansion
# ================================================================================================
    ∇ = f_diff_weights_array(k)
    w = [f_diff_expansion_weights(l[i], ∇) for i ∈ eachindex(l)]
    return append!(repeat(w[1:m],n-k-1),w)
end

# ==============================================================================

@doc raw"""
    lagrangian_interpolation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=1, i=0)

 ``k^{th}``-order lagrangian-interpolation with ``i`` intermediate point of function ``f``
tabulated in normal order at n points, ``f[1],\ \ldots,\ f[n]``.
#### Example:
```
f = [0.0,1,2,3,4,5,6,7]
domain = 0.0..1.0
(X,Y) = lagrangian_interpolation(f, domain; k=2, i=1); println((X,Y))
 (0.0:0.07142857142857142:1.0, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0])
```
"""
function lagrangian_interpolation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=1, i=0)
# ======================================================================================
#   lagrangian (k+1)-point interpolation at i interpolation points
# ======================================================================================
    n = length(f)
    m = i + 1

    l = f_diff_expansion_coeffs_array_interpolation(k, m)
    w1 = f_diff_expansion_weights_array(n, k, m, l)
    w2 = f_diff_function_sequences(f, k, m)
    X = range(domain.left, domain.right, length=(n-1)*m+1)
    Y = w1 .⋅ w2

    return X, Y

end

# ===================================== f_diff_expansion_coeffs_differentiation(k, x) ====

@doc raw"""
    f_diff_expansion_coeffs_differentiation(k::Int, x::T) where T<:Real

Finite-difference expansion coefficient vector ``[l_0(x),\ \ldots\,\ l_p(x)]`` defining
``(k+1)``*-point lagrangian differentiation* of the tabulated analytic function ``f(n+x)``
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
    b = f_diff_expansion_coeffs_interpolation(k, x)

    return polynom_multiplication_coeffs(a, b)[1:k+1]

end

# ================================f_diff_expansion_coeffs_array_differentiation(k, m) ====

@doc raw"""
    f_diff_expansion_coeffs_array_differentiation(k::Int, m::Int)

Finite-difference expansion coefficient vector ``[l_0(x),\ \ldots\,\ l_p(x)]`` defining
``(k+1)``*-point lagrangian differentiation*  of the tabulated analytic function ``f(n+x)``
for ``k*m+1`` x positions ``[-k,\ \cdots,\ 0]`` with ``m-1`` intermediate points.
#### Examples:
```
k = 2; m = 2
o = f_diff_expansion_coeffs_array_differentiation(k,m); println(o)
 [[0.0, 1.0, 0.5], [0.0, 1.0, 0.0], [0.0, 1.0, -0.5], [0.0, 1.0, -1.0], [0.0, 1.0, -1.5]]

m = 1
o = f_diff_expansion_coeffs_array_differentiation(k,m); println(o)
 [[0.0, 1.0, 0.5], [0.0, 1.0, -0.5], [0.0, 1.0, -1.5]]
```
"""
function f_diff_expansion_coeffs_array_differentiation(k::Int, m::Int)

    p = m ≠ 1 ? [i/m-k for i=0:k*m] : [i-k for i=0:k]

    return [f_diff_expansion_coeffs_differentiation(k, p[i]) for i ∈ eachindex(p)]

end

# ================================f_diff_expansion_coeffs_array_differentiation(k, m) ====

@doc raw"""
    lagrangian_differentiation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=1, i=0)

 ``k^{th}``-order lagrangian-differentiation with ``i`` intermediate point of function ``f``
tabulated in normal order at n points, ``f[1],\ \ldots,\ f[n]``.
#### Example:
```
f = [0.0,1,2,3,4,5]
domain = 0.0..5.0
X,Y = lagrangian_differentiation(f, domain; k=2, i = 0); println(X,Y)
 (0.0:1.0:5.0, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
```
"""
function lagrangian_differentiation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=3, i=0)
# ======================================================================================
#   lagrangian (k+1)-point differentiation at i interpolation points
# ======================================================================================
    n = length(f)
    m = i + 1
    l = f_diff_expansion_coeffs_array_differentiation(k, m)
    w1 = f_diff_expansion_weights_array(n, k, m, l)
    w2 = f_diff_function_sequences(f, k, m)
    X = range(domain.left, domain.right, length=(n-1)*m+1)
    Y =  -(n-1)/(domain.right-domain.left) .* (w1 .⋅ w2)

    return X, Y

end
