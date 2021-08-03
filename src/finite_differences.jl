@doc raw"""
    f_diff_weight(k::Int, i::Int)

Weight coefficient
```math
c_{j}^{k}=(-1)^{j}\binom{k}{j},
```
of the ``k^{th}``-order finite difference operator ``\nabla^k`` and corresponding to the function value ``f[n-j]``.
#### Example:
```
k = 5; i = 3
f_diff_weight(k, i)
 -10
```
"""
f_diff_weight(k::Int, i::Int) = iseven(i) ? Base.binomial(k,i) : -Base.binomial(k,i)

# ==============================================================================

@doc raw"""
    f_diff_weights(k::Int)

Weight coefficients ``[c_0^k,\ \ldots,\ c_k^k]`` defining the ``k^{th}``-order finite difference operator,
```math
\nabla^k f[n] = \sum_{j=0}^{k} c_i^kf[n-j],
```
where ``f[n], ...,f[n-k]`` are elements of a tabulated analytic function.
#### Example:
```@docs
k = 3
f_diff_weights(k)
4-element Vector{Int64}:
  1
 -3
  3
 -1
```
"""
f_diff_weights(k::Int) = [f_diff_weight(k, i) for i=0:k]

# ==============================================================================



@doc raw"""
    f_diff_weights_array(kmax::Int)

Collection of weight coefficients ``[c_0^k,\ \ldots,\ c_j^k]`` defining the finite difference operator ``\nabla^j``  ``(0\le j\le k)``.
#### Example:
```@docs
kmax = 3
∇ = f_diff_weights_array(kmax)
4-element Vector{Vector{Int64}}:
 [1]
 [1, -1]
 [1, -2, 1]
 [1, -3, 3, -1]
```
"""
f_diff_weights_array(kmax::Int) = [f_diff_weights(k)  for k=0:kmax]

# ==============================================================================

@doc raw"""
    f_diff_expansion_weights(a, ∇)

Summation weights ``[b_0^k,\ ,\ldots,\ b_k^k]`` corresponding to the expansion coefficients
``[a_0^k,\ ,\ldots,\ a_k^k]`` of a ``k^{th}``-order finite-difference expansion,
```math
\sum_{p=0}^{k}a_{p}\nabla^{p}f[n]=\sum_{j=0}^{k}b_{j}^{k}f[n-j],
```
where ``f[n], ...,f[n-k]`` are elements of a tabulated analytic function.
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
    l = length(coeffs)
    s = [[coeffs[i] * ∇[i][j] for j=1:i] for i=1:l] # s = coeffs .* ∇  but faster
    return [sum([s[i][j] for i=j:l]) for j=1:l]
end

# ==============================================================================

@doc raw"""
    f_diff_expansion_coeffs_interpolation(k::Int, x::T) where T<:Real

Finite-difference expansion coefficients ``l_i(x)`` for lagrangian interpolation of a tabulated
analytic function at offset position ``0\le x\le -k,
```math
f[n+x] =\sum_{p=0}^{k}l_p(x)\nabla^pf[n] = \sum_{j=0}^{k}r_j^k(x)f[n-j],
```
where ``l_0\equiv 0`` and ``l_p(x) = x(x+1)(x+2)\cdots(x+p-1)/p!`` (for ``p=1,\ \ldots,\ k``)
and ``f[n], ...,f[n-k]`` are elements of the tabulated function. The lagrangian interpolation
weights ``[r_0,\ \ldots,\ r_k]`` are calculated with the function `r = f_diff_expansion_weights(l, ∇)`
####
```
k=3
∇ = f_diff_weights_array(k)
x=-1
l = f_diff_expansion_coeffs_interpolation(k,x)
r = f_diff_expansion_weights(l, ∇)
println(l,r)
 [1, -1, 0, 0][0, 1, 0, 0]
```
"""
function f_diff_expansion_coeffs_interpolation(k::Int, x::T) where T<:Real
# ======================================================================================
#   f_difference expansion coefficients for the interpolation interval -k ≤ x ≤ 0
# ======================================================================================
    x > 0 ? error("Error: outside interpolation range (x > 0)") :
    x < -k ? error("Error: outside interpolation range (x < $(-k))") :
    l = ones(T,k+1)
    for i=1:k
        l[i+1] = l[i]*(x+i-1)/i
    end
    return l
end
