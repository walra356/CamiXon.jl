"""
    f_diff_weight(k, j)

Weight coefficient

``c_{j}^{k}=(-1)^{j}\\binom{k}{j},``

of the ``k^{th}``-order finite difference operator ``\nabla^k`` and corresponding to the function value ``f[n-j]``.
#### Example:
```@docs
k = 5; i = 3
f_diff_weight(k, i)
 -10
```
"""
f_diff_weight(k::Int, i::Int) = iseven(i) ? Base.binomial(k,i) : -Base.binomial(k,i)



"""
    f_diff_weights(k)

Weight coefficients ``c_0^k,\\ \\ldots,\\ c_k^k`` corresponding to the function values
``f[n],\\ \\ldots,\ f[n-k]``, for the ``k^{th}``-order finite difference operator ``\nabla^k``.
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



"""
    f_diff_weights_array(kmax)

Weight coefficients ``c_0^k,\\ \\ldots,\\ c_k^k`` corresponding to the function values
``f[n],\\ \\ldots,\ f[n-k]``, for the ``k^{th}``-order finite difference operators
``\\nabla^0,\\ \\ldots,\ \\nabla^k``.
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
