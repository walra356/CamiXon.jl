# =========================== singletons =======================================

struct fwd
end

struct bwd
end

function isforward(notation)

    notation === fwd && return true

    notation === bwd && return false

    error("Error: unknown fdiff notation type")

end

# ==================== fdiff_weight(k, j, notation) ============================

fwd_diff_weight(k::Int, j::Int) = Base.iseven(k+j) ? Base.binomial(k,j) :
                                                    -Base.binomial(k,j)

bwd_diff_weight(k::Int, j::Int) = Base.iseven(j) ? Base.binomial(k,j) :
                                                  -Base.binomial(k,j)

@doc raw"""
    fdiff_weight(k::Int, j::Int, notation=bwd)

Finite difference weight coefficient,

```math
c_{j}^{k}=(-1)^{k+j}\binom{k}{j}.
```
Application:

`fdiff_weight(k,j,fwd)`] `` → c_{k-j}^k``

`fdiff_weight(k,j,bwd)`] `` → c_j^k``

`fdiff_weight(k,j)`] `` → c_j^k``
#### Example:
```
c(k,j,notation) = fdiff_weight(k,j,notation)
k=5
j=3
c(k,j,fwd) == 10
  true

c(k,j,bwd) == -10
  true

c(k,j,fwd) == c(k,k-j)
  true
```
"""
function fdiff_weight(k::Int, j::Int, notation=bwd)

    o = isforward(notation) ? fwd_diff_weight(k, j) : bwd_diff_weight(k, j)

    return o

end





# =============== fdiff_weights(k::Int, notation=fwd) ==========================

@doc raw"""
    fdiff_weights(k::Int, notation=fwd)

Finite difference weights vector defining the ``k^{th}`` *order finite
difference summation weights*.

Application:

[`fdiff_weights(k)`](@ref) `` → \bar{c}^k ≡ [c_k^k,⋯\ c_1^k,\ c_0^k]``

[`fdiff_weights(k,fwd)`](@ref) `` → \bar{c}^k ≡ [c_k^k,⋯\ c_1^k,\ c_0^k]``

[`fdiff_weights(k,bwd)`](@ref) `` → c^k ≡ [c_0^k,\ c_1^k,⋯\ c_k^k]``

where ``c_j^k ← `` [`fdiff_weight(k,j)`](@ref).
#### Example:
```
fdiff_weights0(3,fwd) == [-1, 3, -3, 1]
  true
```
"""
fdiff_weights(k::Int, notation=fwd) = [fdiff_weight(k, j, notation) for j=0:k]

# ==================== fdiff_weights_array(k, j, notation=fwd) =================

@doc raw"""
    fdiff_weights_array(k::Int, notation=fwd)

Finite difference weights vector array defining the *finite
difference summation weights vectors* for the orders 0, 1,⋯ k.

Application:

[`fdiff_weights_array(k)`](@ref) →
`fwd_diffs ≡ ` ``[\bar{c}^0,\ \bar{c}^1,⋯\ \bar{c}^k ]``

[`fdiff_weights_array(k,fwd)`](@ref) →
`fwd_diffs ≡ ` ``[\bar{c}^0,\ \bar{c}^1,⋯\ \bar{c}^k ]``

[`fdiff_weights_array(k,bwd)`](@ref) →
`bwd_diffs ≡ ` ``[c^0,\ c^1,⋯\ c^k ]``

where ``c^k ← `` [`fdiff_weights(k,bwd)`](@ref) and
``\bar{c}^k ← `` [`fdiff_weights(k,fwd)`](@ref)

#### Example:
```
fdiff_weight(3,0,fwd), fdiff_weight(3,0,bwd)
  (-1, 1)

fdiff_weights(3,fwd), fdiff_weights(3,bwd)
  ([-1, 3, -3, 1], [1, -3, 3, -1])

fdiff_weights_array(3,fwd), fdiff_weights_array(3,bwd)
  ([[1], [-1, 1], [1, -2, 1], [-1, 3, -3, 1]], [[1], [1, -1], [1, -2, 1], [1, -3, 3, -1]])
```
"""
fdiff_weights_array(k::Int, notation=fwd) = [fdiff_weights(p, notation)  for p=0:k]

# ============== fdiff_expansion_weights(coeffs, fdiffs, fwd) ==================

@doc raw"""
    fdiff_expansion_weights(coeffs, fdiffs, fwd)

Expansion weights corresponding to the expansion coefficients `coeffs` for
the finite difference expansion `fdiff`.

**Forward difference notation** (`notation = fwd`)

Weight vector ``F^k ≡ [F_k^k,⋯\ F_0^k]`` corresponding to the
expansion coefficients ``α ≡ [α_0^k,⋯\ α_k^k]`` of the ``k^{th}``-order
*forward-difference* expansion,

```math
\sum_{p=0}^{k}α_{p}Δ^{p}f[n]
=\sum_{j=0}^{k}F_{j}^{k}f[n+j]
=F^{k} \cdot f[n:n+k],
```

where ``f[n],⋯\ f[n+k]`` are elements of the
analytic function ``f`` tabulated in *forward* order.

[`fdiff_expansion_weights(coeffs, fdiffs, fwd)`](@ref)
``→ F^k ≡ [F_0^k,⋯\ F_k^k]``,

where `fdiffs ≡ `[`fdiff_weights_array(k)`](@ref) and
`coeffs` = `` α ≡ [α_0,⋯\ α_k]`` defines the expansion.

**Backward difference notation** (`notation = bwd`)

Weight vector ``\bar{B}^{k} ≡ [B_k^k,⋯\ B_0^k]`` corresponding to the
expansion coefficients ``β ≡ [β_0,⋯\ β_k]`` of
the ``k^{th}``-order *backward-difference* expansion,

```math
\sum_{p=0}^{k}β_{p}∇^{p}f[n]
=\sum_{j=0}^{k}B_{k-j}^kf[n-k+j]
=\bar{B}^k \cdot f[n-k:n].
```

where ``f[n-k],⋯\ f[n]`` are elements of the
analytic function ``f`` tabulated in *forward* order.

[`fdiff_expansion_weights(coeffs, fdiffs, bwd)`](@ref)
`` → \bar{B}^{k} ≡ [B_k^k,⋯\ B_0^k]``,

where `fdiffs ≡ `[`fdiff_weights_array(k)`](@ref) and
`coeffs` = `` β ≡ [β_0,⋯\ β_k]`` defines the expansion.
#### Example:
```
k=5
α = β = UnitRange(0,k)
fdiffs = fdiff_weights_array(k)
Fk = fdiff_expansion_weights(α, fdiffs, fwd); println("Fk = $(Fk)")
bBk = fdiff_expansion_weights(β, fdiffs, bwd); println("bBk = $(bBk)")
  Fk = [15, -55, 85, -69, 29, -5]
  bBk = [-5, 29, -69, 85, -55, 15]

bBk == reverse(Fk)
  true
```
"""
function fdiff_expansion_weights(coeffs, fdiffs, notation=fwd)

    forward = isforward(notation)
    c = coeffs
    w = fdiffs
    k = Base.length(c)-1

    o = forward ? [sum([c[1+p] * w[1+p][1+p-j] for p=j:k]) for j=0:k] :
                  [sum([c[1+p] * w[1+p][1+p-j] for p=j:k]) for j=k:-1:0]

end

# ==============================================================================

@doc raw"""
    fdiff_expansion_coeffs_interpolation(k::Int, x::T, fwd) where T<:Real

Finite-difference expansion coefficient vector defining the ``k^{th}``-order
lagrangian interpolation of the tabulated analytic function ``f(n+x)``
at offset position ``x`` with respect to position ``n``.

**Forward difference notation** (`notation = fwd`)

```math
f[n-x] = (1 + ∇)^{-x} f[n] \equiv \sum_{p=0}^k α_p Δ^p f[n] + ⋯,
```
where ``(x)_{p}`` is the [`pochhammer`](@ref) symbol.
Interpolation corresponds to the interval ``-k\le\ x\le 0``;
extrapolation to ``x\ge 0``.

[`fdiff_expansion_coeffs_interpolation(k, x, fwd)`](@ref)
→ ``α^k ≡ [α_0,⋯\ α_k]``

**Backward difference notation** (`notation = bwd`)

```math
f[n+x] = (1 - ∇)^{-x} f[n] \equiv \sum_{p=0}^k β_p ∇^p f[n] + ⋯,
```

[`fdiff_expansion_coeffs_interpolation(k, x, bwd)`](@ref)
→ ``β^k ≡ [β_0,⋯\ β_k]``

#### Examples:
```
k = 5; x = 1
l = fdiff_expansion_coeffs_interpolation(k, x, fwd); println(l)
  [1, -1, 1, -1, 1, -1]

l = fdiff_expansion_coeffs_interpolation(k, x, bwd); println(l)
  [1, 1, 1, 1, 1, 1]
```
"""
function fdiff_expansion_coeffs_interpolation(k::Int, x::T, notation=fwd) where T<:Real

    sgn = notation === fwd ? -1 : notation === bwd ? 1 :
                                  error("Error: unknown notation type")

    l = Base.ones(T,k+1)
    x == 0 ? (for p=2:k+1; l[p] = 0 end) :
                          (for p=1:k; l[p+1] = sgn*l[p]*(x+p-1)/p end)

    return l

end

# ==============================================================================

function lagrange_polynom(f::Vector{T}, x::T, notation=fwd) where T <: Real
# ==============================================================================
#   lagrangian (k+1)-point interpolation at i interpolation points
# ==============================================================================

    Δ = fdiff_weights_array(k)
    α = fdiff_expansion_coeffs_interpolation(k, x, fwd)
    w = fdiff_expansion_weights(α, Δ, fwd)

    return w ⋅ f

end

# =========== fdiff_expansion(coeffs, f, notation=fwd) =================

# ..............................................................................
function fwd_expansion_weights(α)

    k = Base.length(α)-1
    o = [sum([α[p+1] * fdiff_weight(p, j, fwd)  for p=j:k]) for j=0:k]

    return o

end
# ..............................................................................
function bwd_expansion_weights(β)

    k = Base.length(β)-1
    o = [sum([β[p+1] * fdiff_weight(p, j, bwd) for p=j:k]) for j=k:-1:0]

    return o

end
# ------------------------------------------------------------------------------
function fdiff_expansion_weights(coeff, notation=fwd)

    o = isforward(notation) ? fwd_expansion_weights(coeff) :
                              bwd_expansion_weights(coeff)

    return o

end
# ------------------------------------------------------------------------------

@doc raw"""
    fdiff_expansion(coeffs, f, notation=fwd)

Finite difference expansion of the analytical function f(x) tabulated
in *forward order* (growing index) at ``k+1`` positions on a uniform grid.
The expansion coefficients are specified by the vector `coeffs`. By default
`coeffs` are assumed to be in forward-difference notation (`fwd`). For `coeffs`
in backward-difference notation the third argument must be `bwd`.

**Forward difference notation**
```math
\sum_{p=0}^{k}α_{p}Δ^{p}f[n] = F^{k} \cdot f[n:n+k],
```
where ``f[n],⋯\ f[n+k]`` are elements of the
analytic function ``f`` (tabulated in *forward* order) and ``α ≡ [α_0,⋯\ α_k]``
is the vector `coeffs` defining the forward-difference expansion.
The corresponding weights vector ``F^{k}`` is internally generated.

**Backward difference notation**
```math
\sum_{p=0}^{k}β_{p}∇^{p}f[n] = \bar{B}^k \cdot f[n-k:n].
```
where ``f[n-k],⋯\ f[n]`` are elements of the
analytic function ``f`` (tabulated in *forward* order) and
``β ≡ [β_0,⋯\ β_k]`` is the vector `coeffs` defining the backward-difference
expansion. The corresponding weights vector ``\bar{B}^k`` is internally
generated.

#### Examples:
Consider the function ``f(x)=x^2`` and the expansions,
```math
f(x-1)=(1+Δ)^{-1}=(1-Δ+Δ^2-Δ^3+⋯)f(x).
```
```math
f(x+1)=(1-∇)^{-1}=(1+∇+∇^2+∇^3+⋯)f(x),
```
To third order `(k=3)` the forward- and backward-difference coefficient vectors
are `α=[1,-1,1,-1]` and `β=[1,1,1,1]`, respectively. We tabulate the function
at ``k+1`` points, `f=[1,4,9,16]`.
```
α=[1,-1,1,-1]
β=[1,1,1,1]
f=[1,4,9,16]
fdiff_expansion(α, f)      # n=1, f[n]=1, f[n-1] → 0
 0

fdiff_expansion(β, f, bwd) # n=4, f[n]=16, f[n+1] → 25
 25
```
In this case the result is exact because the function is quadratic and
the expansion is third order (lagrangian expansion is based on the polynomial
of ``k^{th}`` degree running through the ``k+1`` points of the tabulated
function).
"""
function fdiff_expansion(coeffs, f, notation=fwd)

    o = fdiff_expansion_weights(coeffs, notation) ⋅ f

    return o

end

# ======================== fdiff_interpolation(f, x; k=3)
@doc raw"""
    fdiff_interpolation(f::Vector{T}, x::V; k=3) where {T <: Real, V <: Real}

Finite difference lagrangian interpolation (by default *third* order) in
between the elements of the analytic function `f` (uniformly tabulated in
*forward* order) for position `x` in *fractional-index units*).
The interpolation points lie on a polynomial curve
(by default *third* degree) running through the tabulated points.
#### Examples:
```
f = [1,2,3,4,5,6,7]
[fdiff_interpolation(f, x; k) for x=1:0.5:7]
  [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]

f = [1,4,9,16,25,36,49]
k=3
[fdiff_interpolation(f, x; k) for x=1:0.5:7]
 [1.0, 2.25, 4.0, 6.25, 9.0, 12.25, 16.0, 20.25, 25.0, 30.25, 36.0, 42.25, 49.0]
```
In this case the result is exact because the function is quadratic and
the expansion is third order (lagrangian expansion is based on the polynomial
of ``k^{th}`` degree running through the ``k+1`` points of the tabulated
function).

![Image](./assets/lagrangian_interpolation.png)
"""
function fdiff_interpolation(f::Vector{T}, x::V; k=3) where {T <: Real, V <: Real}

    l = length(f)
    n = floor(Int,x)
    k = min(k,l-1)
    n = n < l-k ? n : n-k > 1 ? n-k : 1
    α = fdiff_expansion_coeffs_interpolation(k, n-x, fwd)
    o = fdiff_expansion(α, f[n:n+k], fwd)

    return o

end

# ============== fdiff_expansion_coeffs_differentiation(k, x) =================

@doc raw"""
    fdiff_expansion_coeffs_differentiation(k::Int, x::T) where T<:Real

Finite-difference expansion coefficient vector ``[l_0^{\prime}(x),\ ⋯,
\ l_p^{\prime}(x)]`` defining ``k^{th}``-order lagrangian *differentiation*
of the tabulated analytic function ``f(n+x)`` at position ``x``,
```math
\frac{df}{dx}[n+x]=\sum_{p=0}^{k}l_{p}^{\prime}(x)∇^{p}f[n]
```
#### Example:
```
k = 2; x = 0
o = fdiff_expansion_coeffs_differentiation(k,x); println(o)
 [0.0, 1.0, -1.5]
```
"""
function fdiff_expansion_coeffs_differentiation(k::Int, x::T) where T<:Real
# ==============================================================================
#   finite difference expansion coeffs for differentiation
#   in interval -k ≤ x ≤ 0
# ==============================================================================
    a = Base.prepend!([1//i for i=1:k],[0//1])
    b = CamiXon.fdiff_expansion_coeffs_interpolation(k, x, bwd)

    a,b = Base.promote(a,b)

    return CamiXon.polynom_product_expansion(a, b, k)

end

# =================== create_lagrange_differentiation_weights(k, x) ============

@doc raw"""
    create_lagrange_differentiation_weights(k::Int, x::T) where T<:Real

``k^{th}``-order Lagrange differentiation weights vector,
``s^k(x) ≡ [s_k^k(x),⋯\ s_0^k(x)]``, where ``x`` is the position
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

    ∇ = CamiXon.fdiff_weights_array(k)
    coeffs = CamiXon.fdiff_expansion_coeffs_differentiation(k,-k+x)

    return CamiXon.fdiff_expansion_weights(coeffs, ∇, bwd)

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

    ∇ = CamiXon.fdiff_weights_array(k)

    for i=0:k
        coeffs = CamiXon.fdiff_expansion_coeffs_differentiation(k,-k+i)
        m[1+i,1:k+1] = CamiXon.fdiff_expansion_weights(coeffs, ∇, bwd)
    end

    return m

end

# ================ fdiff_expansion_coeffs_array_differentiation(k, m) =========



# =========== trapezoidal_weights(k; rationalize=false, devisor=false) =========

@doc raw"""
    trapezoidal_weights(k::Int [; rationalize=false [, devisor=false]])

Weight coefficient vector ``a=[a_1,⋯\ a_k]`` of trapeziodal rule
optimized for functions of polynomial form,
```math
    ∫_0^n f(x) dx = a_1 (f_0+f_n) + ⋯ + a_k (f_{k-1}+f_{n-k+1})
                                                         + (f_k+⋯+f_{n-k}),
```
where ``k`` is *odd*. The rule is exact for polynonials of degree ``d=0,\ 1,
⋯,\ k-1``. For ``k=1`` the rule reduces to the ordinary trapezoidal rule.
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
  (9, 7257600, [2082753, 11532470, 261166, 16263486, -1020160, 12489922,
                                                     5095890, 7783754, 7200319])
```
"""
function trapezoidal_weights(k::Int; rationalize=false, devisor=false)
# ==============================================================================
# trapezoidal_weights(k; rationalize=false, devisor=false)
# ==============================================================================
    strWarn = "Warning: k = $(k-1) → $(k) (trapezoidal rule requires odd k)"
    Base.isodd(k) ? true :
                   (k=k+1; println(strWarn))

    l = k - 1
    σ = Base.Matrix{Int}(undef,k,k)

    for i=0:k-1
        σ[i+1,:] = CamiXon.polynom_power([i,-1],l) # corresponds to a0,...ak
        σ[i+1,1] = σ[i+1,1] + i^l                  # correction of coeff a0
    end

    F = CamiXon.faulhaber_polynom(k)
    s = [CamiXon.polynom_power([-k,1],p) for p=0:k] .* F
    s[1][1] = -CamiXon.faulhaber_summation(k-1,l)//1

    c = [Base.sum([s[p+1][i+1] for p=i:k]) for i=0:k][1:end-1]

    σ = Base.inv(Base.transpose(σ))

    o = -σ * c  # solving matrix equation results in trapezoidal_weights

    if rationalize
        a = CamiXon.fdiff_expansion_coeffs_adams_moulton(k)
        D = Base.denominator(Base.gcd(a))          # == Adams-Moulton devisor
        o = devisor ? (k, D, Base.round.(Int, o* D)) :
                      Base.round.(Int, o* D) // D  # convert to Rational{Int}
    end

    return o

end

# ======================== trapezoidal_integration(f, domain; k=5) =============

@doc raw"""
    trapezoidal_integration(f, domain, weights)

Integral of the tabulated function ``f=[f_0,⋯\ f_n]`` over the `domain`
``a..b`` using the optimized trapezoidal rule with endpoint correction by the
weightsvector `weights`,
```math
    ∫_0^n f(x) dx = a_1 (f_0+f_n) + ⋯ + a_k (f_{k-1}+f_{n-k+1})
                                                         + (f_k+⋯+f_{n-k}).
```
The rule is exact for polynonials of degree ``d=0,\ 1,⋯\ k-1``.
For ``k=1`` the rule reduces to the ordinary trapezoidal rule (weights = [1/2]).
#### Examples::
```
p = 3
c = [1 for i=0:p]
pol = ImmutablePolynomial(c,:z)
Ipol = integrate(pol)
n = 10

x1=0.0
x2=5.0
x = collect(range(x1, x2, n))
f = pol.(x .-2.5)

w3 = trapezoidal_weights(3)
trapezoidal_integration(f, x1, x2, w3)
 15.416666666666673

Ipol(2.5)-Ipol(-2.5)
 15.41666666666666
```
"""
function trapezoidal_integration(f, x1, x2, weights)

    n = Base.length(f)
    k = Base.length(weights)
    s = (x2-x1)/(n-1)
    a = Base.ones(n); a[1:k] = weights; a[end-k+1:end] = Base.reverse(weights)
    o = (f ⋅ a) * s

    return o

end
