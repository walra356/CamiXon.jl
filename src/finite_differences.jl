@doc raw"""
    f_diff_weight(k, j)

Finite difference weight coefficient
```math
c_{j}^{k}=(-1)^{j}\binom{k}{j},
```

`f_diff_weight(k,j)`] `` \rightarrow c_j^k``
#### Example:
```
c(k, j) = f_diff_weight(k, j)
c(5,3)
 -10
```
"""
f_diff_weight(k::Int, j::Int) = Base.iseven(j) ? Base.binomial(k,j) : -Base.binomial(k,j)

# ==============================================================================

@doc raw"""
    f_diff_weights(k)

Finite difference weights vector ``c^k=[c_k^k,\ \ldots,\ c_0^k]`` defining
the ``k^{th}``-order finite difference operators.

Applications:

**Forward difference notation**

The *forward difference* summation is
```math
\Delta^k f[n]=[c_{k}^{k},\thinspace c_{k-1}^{k},\thinspace\ldots,c_{0}^{k}]\left[\begin{array}{c}
f[n]\\
\vdots\\
f[n+k]
\end{array}\right]=\sum_{j=0}^{k} c_{k-j}^kf[n+j].
```

This form is designed for use with *analytical* functions, ``f``, tabulated
in *forward* order as ``f[n], ...,f[n+k]``.

**Backward difference notation**

The *backward difference* summation is
```math
\nabla^{k}f[n]=[c_{k}^{k},\thinspace c_{k-1}^{k},\thinspace\ldots,c_{0}^{k}]\left[\begin{array}{c}
f[n-k]\\
\vdots\\
f[n]
\end{array}\right]=\sum_{j=0}^{k}c_{k-j}^kf[n-k+j].
```

This form is designed for use with *analytical* functions, ``f``, tabulated
in *forward* order as ``f[n-k], ...,f[n]``.

`f_diff_weights(k)` `` \rightarrow \ c^k ≡ [c_k^k,\ c_1^k,\ldots,\ c_0^k]``,

where [`f_diff_weight(k,j)`](@ref) `` \rightarrow c_j^k``.
#### Example:
```
c(k) = f_diff_weights(k)
c(3)
4-element Vector{Int64}:
  1
 -3
  3
 -1
```
"""
f_diff_weights(k::Int) = [CamiXon.f_diff_weight(k, k-j) for j=0:k]

# ==============================================================================

@doc raw"""
    f_diff_weights_array(kmax)

Collection of finite difference weight vectors, ``c^0,\ \ldots,\ c^k``, where
``c^k`` = [`f_diff_weights(k)`](@ref).

Application in [`Finite difference expansions`](@ref).

`f_diff_weights_array(kmax)` ``\rightarrow\ [\ c^0,\ c^1,\ \ldots,\ c^{kmax} ]``,

where [`f_diff_weights(k)`](@ref)``\rightarrow\ c^k ≡ [c_k^k,\ c_1^k,\ldots,\ c_0^k]``.
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
f_diff_weights_array(kmax::Int) = [CamiXon.f_diff_weights(k)  for k=0:kmax]

# ==============================================================================

@doc raw"""
    f_diff_expansion_weights(α, ∇)

Weight vector ``b^k ≡ [b_k^k,\ ,\ldots,\ b_0^k]`` corresponding to the
expansion coefficients ``[α_0^k,\ ,\ldots,\ α_k^k]`` of the ``k^{th}``-order
finite-difference expansion,

```math
\sum_{p=0}^{k}α_{p}\nabla^{p}f[n]=\sum_{j=0}^{k}b^k[j]f[n-k+j],
```

where ``b^k[j] \equiv b_{k-j}^k`` and ``f[n-k], ...,f[n]`` are elements of the
analytic function ``f`` tabulated in *forward* order. Note the difference in
ordering between the finite-difference expansion *coefficients*,
``α_{0},\ \ldots,\ α_{k}``, and the finite-difference expansion *weights*,
``b_k^{k},\ \ldots,\ b_0^{k}``. Note further the difference in ``k`` dependence:
the *weights*, ``b_j^k``, are ``k``*-dependent*, whereas the *coefficients*,
``α_j``, are not.
#### Example:
```
k=5
∇ = f_diff_weights_array(k)
α = UnitRange(0,k)
b = f_diff_expansion_weights(α, ∇)
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
    k = Base.length(coeffs)-1
    #return [sum([coeffs[1+p] * ∇[1+p][1+p-j] for p=j:k]) for j=0:k]
    return [Base.sum([coeffs[1+p] * ∇[1+p][1+p-j] for p=j:k]) for j=k:-1:0]
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
    l = Base.ones(T,k+1)
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
``f`` tabulated in forward order on a uniform grid of ``n`` points, ``f[1], ...,f[n]``;
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
    n = Base.length(f)
    return [f[CamiXon.summation_range(n,i,k,m)] for i=0:(n-1)*m]
end

# ==============================================================================

@doc raw"""
    lagrangian_interpolation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=1, m=1)

``k^{th}``-order lagrangian *interpolation* of the analytic function ``f``
tabulated in forward order on a uniform grid of ``n`` points, ``f[1],\ \ldots,
\ f[n]``; ``m`` is the multiplier defining the interpolation grid size.
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

    ∇ = CamiXon.f_diff_weights_array(k)
    l = [CamiXon.f_diff_expansion_coeffs_lagrange(k, x) for x=-k:1/m:0]
    w = [CamiXon.f_diff_expansion_weights(l[i], ∇) for i ∈ eachindex(l)]
    w1 = Base.append!(Base.repeat(w[1:m],n-k-1),w)
    w2 = CamiXon.f_diff_function_sequences(f, k, m)

    X = Base.range(domain.left, domain.right, length=(n-1)*m+1)
    Y = [w1[i] ⋅ w2[i] for i ∈ Base.eachindex(w1)]

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
    n = Base.length(f)

    ∇ = CamiXon.f_diff_weights_array(k)
    l = [CamiXon.f_diff_expansion_coeffs_lagrange(k, x) for x=0:1/m:e]
    w1 = [CamiXon.f_diff_expansion_weights(l[i], ∇) for i ∈ Base.eachindex(l)]
    w2 = CamiXon.f_diff_function_sequences(f, k, m)[end]

    ΔX = (domain.right - domain.left)/((n-1)*m)
    X = Base.range(domain.right, domain.right + ΔX * m*e, length=m*e+1)
    Y = [w1[i] ⋅ w2 for i ∈ Base.eachindex(w1)]

    return X, Y

end

# ===================================== f_diff_expansion_coeffs_differentiation(k, x) ====

@doc raw"""
    f_diff_expansion_coeffs_differentiation(k::Int, x::T) where T<:Real

Finite-difference expansion coefficient vector ``[l_0^{\prime}(x),\ \ldots,
\ l_p^{\prime}(x)]`` defining ``k^{th}``-order lagrangian *differentiation*
of the tabulated analytic function ``f(n+x)`` at position ``x``,
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
    a = Base.prepend!([1//i for i=1:k],[0//1])
    b = CamiXon.f_diff_expansion_coeffs_lagrange(k, x)

    a,b = Base.promote(a,b)

    return CamiXon.polynom_product_expansion(a, b, k)

end

# =================== create_lagrange_differentiation_weights(k, x) ============

@doc raw"""
    create_lagrange_differentiation_weights(k::Int, x::T) where T<:Real

``k^{th}``-order Lagrange differentiation weights vector,
``s^k(x) ≡ [s_k^k(x),\ ,\ldots,\ s_0^k(x)]``, where ``x`` is the position
relative point ``n``.

```math
\frac{df}{dx}[n+x]= \sum_{j=0}^{k}s_{k-j}^k(x)f[n-k+j],
```
where ``s^k_x[j] ≡ s_{k-j}(x)^k``.
#### Example:
```
k = 3
x = 0
ldw = create_lagrange_differentiation_weights(k,x); println(ldw)
  Rational{Int64}[-11//6, 3//1, -3//2, 1//3]

 sum(ldw)
   0//1
```
"""
function create_lagrange_differentiation_weights(k::Int, x::T) where T<:Real

    ∇ = CamiXon.f_diff_weights_array(k)
    coeffs = CamiXon.f_diff_expansion_coeffs_differentiation(k,-k+x)

    return CamiXon.f_diff_expansion_weights(coeffs,∇)

end

# ============ create_lagrange_differentiation_matrix(k) =======================

@doc raw"""
    create_lagrange_differentiation_matrix(k::Int)

Lagrange differentiation matrix, ``m[i,j]=s_{k-j}^k(i)``, for ``k^{th}``-order
lagrangian differentiation,
```math
\frac{dy}{dx}[i]= \sum_{j=0}^{k}m[i,j]y[j],
```
#### Example:
```
k = 3
create_lagrange_differentiation_matrix(k)
 4×4 Matrix{Rational{Int64}}:
  -11//6   3//1  -3//2   1//3
   -1//3  -1//2   1//1  -1//6
    1//6  -1//1   1//2   1//3
   -1//3   3//2  -3//1  11//6
```
"""
function create_lagrange_differentiation_matrix(k::Int)

    m = Matrix{Rational{Int}}(undef,k+1,k+1)

    ∇ = CamiXon.f_diff_weights_array(k)

    for i=0:k
        coeffs = CamiXon.f_diff_expansion_coeffs_differentiation(k,-k+i)
        m[1+i,1:k+1] = CamiXon.f_diff_expansion_weights(coeffs,∇)
    end

    return m

end

# ================ f_diff_expansion_coeffs_array_differentiation(k, m) =========

@doc raw"""
    lagrange_differentiation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=1, m=1)

``k^{th}``-order lagrangian *differentiation* of the analytic function ``f``,
tabulated in forward order on a uniform grid of ``n`` points, ``f[1],\ \ldots,
\ f[n]``; ``m`` is the multiplier for intermediate positions (for ``m=1``
*without* intermediate points).
#### Example:
```
f = [0.0,1,4,9,16,25] # f = x^2
domain = 0.0..5.0
X,Y = lagrange_differentiation(f, domain; k=2, m = 1)
  (0.0:1.0:5.0, [0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
```
"""
function lagrange_differentiation(f::Vector{Float64}, domain::ClosedInterval{Float64}; k=1, m=1)
# ==============================================================================
#   lagrangian (k+1)-point differentiation at m interpolation points
# ==============================================================================
    n = Base.length(f)

    ∇ = CamiXon.f_diff_weights_array(k)
    l = [CamiXon.f_diff_expansion_coeffs_differentiation(k, x) for x=-k:1/m:0]
    w = [CamiXon.f_diff_expansion_weights(l[i], ∇) for i ∈ Base.eachindex(l)]
    w1 = Base.append!(repeat(w[1:m],n-k-1),w)
    w2 = CamiXon.f_diff_function_sequences(f, k, m)

    X = Base.range(domain.left, domain.right, length=(n-1)*m+1)
    Y = (n-1)/(domain.right-domain.left) .*  [w1[i] ⋅ w2[i] for i ∈ Base.eachindex(w1)]

    return X, Y

end

# =========== trapezoidal_weights(k; rationalize=false, devisor=false) =========

@doc raw"""
    trapezoidal_weights(k::Int [; rationalize=false [, devisor=false]])

Weight coefficient vector ``a=[a_1,\cdots,\ a_k]`` of trapeziodal rule
optimized for functions of polynomial form,
```math
    ∫_0^n f(x) dx = a_1 (f_0+f_n)+\cdots+a_k (f_{k-1}+f_{n-k+1}) + (f_k+\cdots+f_{n-k}),
```
where ``k`` is *odd*. The rule is exact for polynonials of degree ``d=0,\ 1,
\cdots,\ k-1``. For ``k=1`` the rule reduces to the ordinary trapezoidal rule.
By default the output is in Float64, optionally the output is rational, with or
without specification of the gcd devisor.
#### Example::
```
[trapezoidal_weights(k; rationalize=true, devisor=true) for k=1:2:9]
5-element Vector{Tuple{Int64, Int64, Vector{Int64}}}:
  (1, 2, [1])
  (3, 24, [9, 28, 23])
  (5, 1440, [475, 1902, 1104, 1586, 1413])
  (7, 120960, [36799, 176648, 54851, 177984, 89437, 130936, 119585])
  (9, 7257600, [2082753, 11532470, 261166, 16263486, -1020160, 12489922, 5095890, 7783754, 7200319])
```
"""
function trapezoidal_weights(k::Int; rationalize=false, devisor=false)
# ==============================================================================
# trapezoidal_weights(k; rationalize=false, devisor=false)
# ==============================================================================
    Base.isodd(k) ? true : (k=k+1; println("Warning: k = $(k-1) → $(k) (trapezoidal rule requires odd k)"))

    l = k - 1
    σ = Base.Matrix{Int}(undef,k,k)

    for i=0:k-1
        σ[i+1,:] = CamiXon.polynom_power([i,-1],l)      # corresponds to coeffs a0,...,ak
        σ[i+1,1] = σ[i+1,1] + i^l                       # correction of coeff a0
    end

    F = CamiXon.faulhaber_polynom(k)
    s = [CamiXon.polynom_power([-k,1],p) for p=0:k] .* F
    s[1][1] = -CamiXon.faulhaber_summation(k-1,l)//1

    c = [Base.sum([s[p+1][i+1] for p=i:k]) for i=0:k][1:end-1]

    σ = Base.inv(Base.transpose(σ))

    o = -σ * c  # solving matrix equation results in trapezoidal_weights as real numbers

    if rationalize
        a = CamiXon.f_diff_expansion_coeffs_adams_moulton(k)
        D = Base.denominator(Base.gcd(a))       # trapezoidal_weights divisor == Adams-Moulton devisor
        o = devisor ? (k, D, Base.round.(Int, o* D)) : Base.round.(Int, o* D) // D              # convert to Rational{Int}
    end

    return o

end

# ======================== trapezoidal_integration(f, domain; k=5) =============

@doc raw"""
    trapezoidal_integration(f, domain, weights)

Integral of the tabulated function ``f=[f_0,\cdots,\ f_n]`` over the `domain`
``a..b`` using the optimized trapezoidal rule with endpoint correction by the
weightsvector `weights`,
```math
    ∫_0^n f(x) dx = a_1 (f_0+f_n)+\cdots+a_k (f_{k-1}+f_{n-k+1}) + (f_k+\cdots+f_{n-k}).
```
The rule is exact for polynonials of degree ``d=0,\ 1,\cdots,\ k-1``.
For ``k=1`` the rule reduces to the ordinary trapezoidal rule (weights = [1/2]).
#### Examples::
```
p = 3
c = [1 for i=0:p]
pol = ImmutablePolynomial(c,:z)
Ipol = integrate(pol)
n = 10

domain = 0.0..5.0
x = collect(range(domain, n))
f = pol.(x .-2.5)

w3 = trapezoidal_weights(3)
trapezoidal_integration(f, domain, w3)
 15.416666666666673

Ipol(2.5)-Ipol(-2.5)
 15.41666666666666
```
"""
function trapezoidal_integration(f, domain, weights)

    n = Base.length(f)
    k = Base.length(weights)
    s = (domain.right-domain.left)/(n-1)
    a = Base.ones(n); a[1:k] = weights; a[end-k+1:end] = Base.reverse(weights)
    o = (f ⋅ a) * s

    return o

end
