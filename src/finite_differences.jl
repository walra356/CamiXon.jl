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



@doc raw"""
    f_diff_weights(k::Int)

Weight coefficients ``[c_0^k,\ \ldots,\ c_k^k]`` defining the ``k^{th}``-order finite difference operator,
```math
\nabla^k f[n] = \sum_{j=0}^{k} c_i^kf[n-j].
```
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



@doc raw"""
    f_diff_weights_array(kmax::Int)

Collection of weight coefficients ``c_0^k,\ \ldots,\ c_k^j`` defining the finite difference operators ``\nabla^j`` with ``j\in0\ \ldots,\ k``.
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

@doc raw"""
    f_diff_expansion_weights(a, ∇)

Summation weights ``[b_0^k,\ ,\ldots,\ b_k^k]`` corresponding to the expansion coefficients
``[a_0^k,\ ,\ldots,\ a_k^k]`` of a ``k^{th}``-order finite-difference expansion,

```math
\sum_{p=0}^{k}a_{p}\nabla^{p}f[n]=\sum_{j=0}^{k}b_{j}^{k}f[n-j],
```

where ``f[n], ...,f[n-k]`` are ``k+1`` values of a tabulated anaytic function f.
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
