@doc raw"""
    f_diff_weight(k, j)

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
f_diff_weights(k::Int) = [f_diff_weight(k, k-i) for i=0:k]


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

where ``f[n-k], ...,f[n]`` are elements of the analytic function ``f`` tabulated in *normal ordering*.
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
####
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
    summation_range(n::Int, j::Int, k::Int, i::Int)

*Finite-difference summation range* for position ``j\in (0,\ 1,\ \ldots,\ n(i+1))`` as used
in ``(k+1)``*-point lagrangian interpolation* with ``i`` intermediate points.
#### Examples:
```
n = 7; j = 2, k = 3; i = 0
a = summation_range(n, j, k, i); println(a)
b = [summation_range(n, j, k, i) for j=0:(n-1)*m]; println(b)
 3:6
 UnitRange{Int64}[1:4, 2:5, 3:6, 4:7, 4:7, 4:7, 4:7]
```
"""
function summation_range(n::Int, j::Int, k::Int, i::Int) 
# ================================================================================================
#   summation range for position j in lagrangian interpolation with i intermediate points
# ================================================================================================
    m = i + 1
    return j < (n-1-k)*m  ? UnitRange(j÷m+1,j÷m+k+1) : UnitRange(n-k,n)
end
