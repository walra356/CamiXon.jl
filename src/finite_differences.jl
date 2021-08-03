"""
    f_diff_weight(k, j)

Weight coefficient

``c_{j}^{k}=(-1)^{j}\binom{k}{j},``

of the ``k^{th}``-order finite difference operator ``\nabla^k`` and corresponding to the function value ``f[n-j]``.
#### Example:
```@docs
k = 5; i = 3
f_diff_weight(k, i)
 -10
```
"""
f_diff_weight(k::Int, i::Int) = iseven(i) ? Base.binomial(k,i) : -Base.binomial(k,i)
