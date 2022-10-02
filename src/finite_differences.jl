# ==================== fdiff_weight(k, j) ======================================
@doc raw"""
    fdiff_weight(k::Int, j::Int)

Finite difference weight coefficient,

```math
c_{j}^{k}=(-1)^{k+j}\binom{k}{j}.
```
#### Example:
```
c(k,j) = fdiff_weight(k,j)

[[c(k,j) for j=0:k] for k=0:3] == [[1], [1, -1], [1, -2, 1], [1, -3, 3, -1]]
  true

[[c(k,k-j) for j=0:k] for k=0:3] == [[1], [-1, 1], [1, -2, 1], [-1, 3, -3, 1]]
  true
```
"""
function fdiff_weight(k::Int, j::Int)

    return Base.iseven(j) ? Base.binomial(k,j) : -Base.binomial(k,j)

end

# ============== fdiff_expansion_weights(coeffs, bwd, rev) =====================

# ..............................................................................
function reg_fwd_expansion_weights(coeffs)

    k = Base.length(coeffs)-1

    o = [Base.sum([coeffs[p+1] * fdiff_weight(p, p-j) for p=j:k]) for j=0:k]

    return o

end
# ..............................................................................
function rev_fwd_expansion_weights(coeffs)

    k = Base.length(coeffs)-1

    o = [Base.sum([coeffs[p+1] * fdiff_weight(p, p-j) for p=j:k]) for j=k:-1:0]

    return o

end
# ..............................................................................
function reg_bwd_expansion_weights(coeffs)

    k = Base.length(coeffs)-1

    o = [sum([coeffs[p+1] * fdiff_weight(p, j) for p=j:k]) for j=0:k]

    return o

end
# ..............................................................................
function rev_bwd_expansion_weights(coeffs)

    k = Base.length(coeffs)-1

    o = [sum([coeffs[p+1] * fdiff_weight(p, j) for p=j:k]) for j=k:-1:0]

    return o

end
# ..............................................................................
@doc raw"""
    fdiff_expansion_weights(coeffs[, notation=bwd[, ordering=rev]])

Expansion weights corresponding to the expansion coefficients `coeffs` of
a finite difference expansion.

**Forward difference notation** (`notation = fwd`)

Weight vector ``F^k ≡ [F_k^k,⋯\ F_0^k]`` corresponding to the
expansion coefficients ``α ≡ [α_0^k,⋯\ α_k^k]`` of the ``k^{th}``-order
*forward-difference* expansion,

```math
\sum_{p=0}^{k}α_{p}Δ^{p}f[n]
=\sum_{j=0}^{k}F_{j}^{k}f[n+j]
=F^{k} \cdot f[n:n+k],
```

where ``f[n:n+k]`` are elements of the
analytic function ``f`` tabulated in *forward* order.

[`fdiff_expansion_weights(α, fwd, reg)`](@ref)
``→ F^k ≡ [F_0^k,⋯\ F_k^k]``,

where `` α ≡ [α_0,⋯\ α_k]`` has to be supplied in combination with `fwd` to
indicate that the weights must be evaluated in forward-difference notation.

**Backward difference notation** (`notation = bwd`)

Weight vector ``\bar{B}^{k} ≡ [B_k^k,⋯\ B_0^k]`` corresponding to the
expansion coefficients ``β ≡ [β_0,⋯\ β_k]`` of
the ``k^{th}``-order *backward-difference* expansion,

```math
\sum_{p=0}^{k}β_{p}∇^{p}f[n]
=\sum_{j=0}^{k}B_{j}^kf[n-j]
=\bar{B}^k \cdot f[n-k:n].
```

where ``f[n-k:n]`` are elements of the
analytic function ``f`` tabulated in *forward* order.

[`fdiff_expansion_weights(β, bwd, rev)`](@ref)
`` → \bar{B}^{k} ≡ [B_k^k,⋯\ B_0^k]``,

where `` β ≡ [β_0,⋯\ β_k]`` has to be supplied in combination with `bwd` to
indicate that the weights must be evaluated in backward-difference notation.
#### Example:

Consider the expansions,

```math
f[n-1]=(1+Δ)^{-1}f[n]=(1-Δ+Δ^2-Δ^3+⋯)f[n].
```
```math
f[n+1]=(1-∇)^{-1}f[n]=(1+∇+∇^2+∇^3+⋯)f[n],
```
```
α = [1,-1,1,-1,1]
β = [1,1,1,1,1]
Fk = fdiff_expansion_weights(α, fwd, reg); println("Fk = $(Fk)")
  Fk = [5, -10, 10, -5, 1]

Bk = fdiff_expansion_weights(β, bwd, reg); println("Bk = $(Bk)")
  Bk = [5, -10, 10, -5, 1]

revFk = fdiff_expansion_weights(α, fwd, rev); println("revFk = $(revFk)")
  revFk = [1, -5, 10, -10, 5]

revBk = fdiff_expansion_weights(β, bwd, rev); println("revBk = $(revBk)")
  revBk = [1, -5, 10, -10, 5]
```
"""
function fdiff_expansion_weights(coeffs, notation=bwd, ordering=rev)

    if isforward(notation)

        o = isregular(ordering) ? reg_fwd_expansion_weights(coeffs) :
                                  rev_fwd_expansion_weights(coeffs)
    else

        o = isregular(ordering) ? reg_bwd_expansion_weights(coeffs) :
                                  rev_bwd_expansion_weights(coeffs)
    end

    return o

end

# =========== fdiff_expansion(coeffs, f, notation=bwd) =========================

# ..............................................................................
@doc raw"""
    fdiff_expansion(coeffs, f[, notation=bwd])

Finite difference expansion of the analytical function f(x) tabulated
in *forward order* (growing index) at ``k+1`` positions on a uniform grid.
The expansion coefficients are specified by the vector `coeffs`. By default
`coeffs` are assumed to be in backward-difference notation (`bwd`). For `coeffs`
in forward-difference notation the third argument must be `fwd`.

**Forward difference notation** (`notation = fwd`)
```math
\sum_{p=0}^{k}α_{p}Δ^{p}f[n] = F^{k} \cdot f[n:n+k],
```
where ``f[n:n+k]`` are elements of the
analytical function ``f`` (tabulated in *forward* order) and
``α ≡ [α_0,⋯\ α_k]`` is the vector `coeffs`, which has to be supplied to
define the forward-difference expansion.
The corresponding weights vector ``F^{k}`` is internally generated.

**Backward difference notation** (`notation = bwd`)
```math
\sum_{p=0}^{k}β_{p}∇^{p}f[n] = \bar{B}^k \cdot f[n-k:n].
```
where ``f[n-k:n]`` are elements of the
analytical function ``f`` (tabulated in *forward* order) and
``β ≡ [β_0,⋯\ β_k]`` is the vector `coeffs`, which has to be supplied to
define the backward-difference expansion. The corresponding weights vector
``\bar{B}^k`` is internally generated.

#### Examples:
Consider the function ``f(x)=x^2`` and the expansions,
```math
f(x-1)=(1+Δ)^{-1}=(1-Δ+Δ^2-Δ^3+⋯)f(x).
```
```math
f(x+1)=(1-∇)^{-1}=(1+∇+∇^2+∇^3+⋯)f(x),
```
To fourth order `(k=4)` the forward- and backward-difference coefficient vectors
are `α=[1,-1,1,-1,1]` and `β=[1,1,1,1,1]`, respectively. We tabulate the function
at ``k+1`` points, `f=[1,4,9,16,25]`.
```
α = [1,-1,1,-1,1]
Fk = fdiff_expansion_weights(α, fwd, reg); println("Fk = $(Fk)")
  Fk = [5, -10, 10, -5, 1]

β = [1,1,1,1,1]
revBk = fdiff_expansion_weights(β, bwd, rev); println("revBk = $(revBk)")
  revBk = [1, -5, 10, -10, 5]

f = [1,4,9,16,25]
o = fdiff_expansion(α, f, fwd); println("f[0] = $(o)")
  f[0] = 0

fdiff_expansion(α, f, fwd) == Fk ⋅ f == fdiff_interpolation(f, 0)
  true

o = fdiff_expansion(β, f, bwd); println("f[6] = $(o)")
  f[6] = 36

fdiff_expansion(β, f, bwd) == revBk ⋅ f == fdiff_interpolation(f, length(f)+1)
  true
```
In these cases the results are exact because the function is quadratic and
the expansion is third order (based on the polynomial of ``k^{th}`` degree
running through the ``k+1`` points of the tabulated function). Note the
relation with [`fdiff_interpolation(f, v, k=3)`](@ref).
"""
function fdiff_expansion(coeffs, f, notation=bwd)

    ordering = isforward(notation) ? reg : rev
    w = fdiff_expansion_weights(coeffs, notation, ordering)

    return w ⋅ f

end

# ======= fdiff_interpolation_expansion_coeffs(k, x, notation=bwd) =============

#...............................................................................
function fwd_interpolation_expansion_coeffs(ξ::T, k=3) where T<:Real

    o = Base.ones(T,k+1)
    ξ == 0 ? (for p=2:k+1; o[p] = 0 end) :
             (for p=1:k; o[p+1] = o[p]*(ξ-p+1)/p end)

    return o

end
#...............................................................................
function bwd_interpolation_expansion_coeffs(ξ::T, k=3) where T<:Real

    o = Base.ones(T,k+1)
    ξ == 0 ? (for p=2:k+1; o[p] = 0 end) :
             (for p=1:k; o[p+1] = o[p]*(p-ξ-1)/p end)

    return o

end
#...............................................................................
@doc raw"""
    fdiff_interpolation_expansion_coeffs(ξ::T [, k=3 [, notation=bwd]]) where T<:Real

Finite-difference expansion coefficient vector defining the ``k^{th}``-order
(default *third* order) Lagrange-polynomial interpolation of a tabulated
analytic function ``f[n]`` at offset ``ξ`` with respect to index
position ``n``, which is positive for increasing index and negative for
decreasing index.

**Forward difference notation** (`notation = fwd`)

In this case we consider the tabulated interval ``f[n:n+k]``. The interpolated
value ``f[n+ξ]`` is given by the forward-difference expansion

```math
f[n+ξ] = \sum_{p=0}^k α_p(-ξ) Δ^p f[n] + ⋯,
```
where the expansion coefficients are given by

[`fdiff_interpolation_expansion_coeffs(ξ, k, fwd)`](@ref)
`` → α(-ξ) ≡ [α_0(-ξ),⋯\ α_k(-ξ)]``. In this notation the range
``0\leq ξ\leq k`` corresponds to interpolation and the ranges ``ξ<0`` and
``ξ>k`` to extrapolation.

**Backward difference notation** (`notation = bwd`)

In this case we consider the tabulated interval ``f[n-k:n]``. The interpolated
value ``f[n+ξ]`` is given by the backward-difference expansion

```math
f[n+ξ] = \sum_{p=0}^k β_p(ξ) ∇^p f[n] + ⋯,
```
where the expansion coefficients are given by

[`fdiff_interpolation_expansion_coeffs(ξ, k, bwd)`](@ref)
`` → β(ξ) ≡ [β_0(ξ),⋯\ β_k(ξ)]``. In this notation the range
``-k\leq ξ\leq0`` corresponds to interpolation and the ranges
``ξ<-k`` and ``ξ>0`` to extrapolation.

#### Examples:
```
k = 5
ξ = -1
α = fdiff_interpolation_expansion_coeffs(ξ, k, fwd); println("α = $α")
β = fdiff_interpolation_expansion_coeffs(ξ, k, bwd); println("β = $β")
  α = [1, 1, 0, 0, 0, 0]
  β = [1, 1, 1, 1, 1, 1]

ξ = 0
α = fdiff_interpolation_expansion_coeffs(ξ, k, fwd); println("α = $α")
β = fdiff_interpolation_expansion_coeffs(ξ, k, bwd); println("β = $β")
  α = [1, 0, 0, 0, 0, 0]
  β = [1, 0, 0, 0, 0, 0]

ξ = 1
α = fdiff_interpolation_expansion_coeffs(ξ, k, fwd); println("α = $α")
β = fdiff_interpolation_expansion_coeffs(ξ, k, bwd); println("β = $β")
  α = [1, -1, 1, -1, 1, -1]
  β = [1, -1, 0, 0, 0, 0]
```
"""
function fdiff_interpolation_expansion_coeffs(ξ::T, k=3, notation=bwd) where T<:Real

    o = isforward(notation) ? fwd_interpolation_expansion_coeffs(-ξ, k) :
                              bwd_interpolation_expansion_coeffs(ξ, k)
    return o

end

# ============== fdiff_interpolation_expansion_weights(coeffs) =================

#...............................................................................
function fwd_interpolation_expansion_weights(ξ::T, k=3, ordering=reg) where T<:Real

    α = fdiff_interpolation_expansion_coeffs(ξ, k, fwd)
    o = fdiff_expansion_weights(α, fwd, ordering)

    return o

end
#...............................................................................
function bwd_interpolation_expansion_weights(ξ::T, k=3, ordering=rev) where T<:Real

    β = fdiff_interpolation_expansion_coeffs(ξ, k, bwd)
    o = fdiff_expansion_weights(β, bwd, ordering)

    return o

end
#...............................................................................
function fdiff_interpolation_expansion_weights(ξ::T, k=3, notation=bwd, ordering=rev) where T<:Real

    o = isforward(notation) ? fwd_interpolation_expansion_weights(-ξ, k, ordering) :
                              bwd_interpolation_expansion_weights(ξ, k, ordering)
    return o

end
function fdiff_interpolation_expansion_weights(coeffs, notation=bwd, ordering=rev)

    if isforward(notation)

        o = isregular(ordering) ? reg_fwd_expansion_weights(coeffs) :
                                  rev_fwd_expansion_weights(coeffs)
    else

        o = isregular(ordering) ? reg_bwd_expansion_weights(coeffs) :
                                  rev_bwd_expansion_weights(coeffs)
    end

    return o

end

# ======================== fdiff_interpolation(f, x; k=3) ======================
@doc raw"""
    fdiff_interpolation(f::Vector{T}, v::V; k=3) where {T<:Real, V<:Real}

Finite difference lagrangian interpolation (by default *third* order) at real
position v (in index units) with respect to the elements of the uniformly
tabulated analytic function `f[1:N]`. The interpolation points lie on a
Lagrange polynomial of degree ``k`` (by default *third* degree) running through
``k+1`` subsequenct points of the tabulated function. Outside the tabulated
range, the method represents extrapolation on the lagrangian polynomial defined
by the first/last ``k+1`` tabulated points.

Beware that the interpolation becomes inaccurate if the tabulated function
cannot be approximated by a polynomial of degree ``k``.
#### Examples:
```
f = [1,2,3,4,5,6,7]
[fdiff_interpolation(f, v; k=3) for v=1:0.5:7]
  [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]

f = [1,4,9,16,25,36,49]
[fdiff_interpolation(f, v; k=3) for v=1:0.5:7]
 [1.0, 2.25, 4.0, 6.25, 9.0, 12.25, 16.0, 20.25, 25.0, 30.25, 36.0, 42.25, 49.0]

 f = [x^3 for x=-4:2]
 f1(v) = fdiff_interpolation(f, v; k=1)
 f2(v) = fdiff_interpolation(f, v; k=2)
 f3(v) = fdiff_interpolation(f, v; k=3)
 [[f1(v),f2(v),f3(v)] for v=1:0.5:9]
   [[-64.0, -64.0, -64.0], [-45.5, -43.25, -42.875], [-27.0, -27.0, -27.0],
   [-17.5, -16.0, -15.625], [-8.0, -8.0, -8.0], [-4.5, -3.75, -3.375],
   [-1.0, -1.0, -1.0], [-0.5, -0.5, -0.125], [0.0, 0.0, 0.0],
   [0.5, -0.25, 0.125], [1.0, 1.0, 1.0], [4.5, 3.75, 3.375], [8.0, 8.0, 8.0],
   [11.5, 13.75, 15.625], [15.0, 21.0, 27.0], [18.5, 29.75, 42.875],
   [22.0, 40.0, 64.0]]
```
The result for f3(v) is exact because the function is cubic and
the expansion is third order - see Figure below. The tabulated function is
given by the black points. The interpolation and extrapolation points are red.

![Image](./assets/lagrangian_interpolation.png)
"""
function fdiff_interpolation(f::Vector{T}, v::V; k=3) where {T<:Real, V<:Real}

    l = length(f)
    k = min(k,l-1)
    k > 0 || error("Error: k ≥ 1 required for lagrangian interpolation")
    n = v < 1 ? 1 : v < l-k ? floor(Int,v) : l-k
    α = fdiff_interpolation_expansion_coeffs(n-v, k, fwd)
    o = fdiff_expansion(α, f[n:n+k], fwd)

    return o

end

# ============== fdiff_differentiation_expansion_coeffs(ξ, k) =================

#...............................................................................
function fwd_differentiation_expansion_coeffs(ξ::T, k=3) where T<:Real
# ==============================================================================
#   Forward difference expansion coefficients for differentiation of an
#   analytic function f(x) tabulated under the convention f[n,n+k] and
#   evaluated at the interpolation position n+ξ.
# ==============================================================================
    Float = (Float64,Float32,Float16)

    a = T ∈ Float ? Base.ones(T,k+1) :
                    T <: Rational{}  ? Base.ones(T,k+1) :
                                       Base.ones(Rational{T},k+1)
    for i=1:k
        a[i+1] = iseven(i) ? -i : i
    end

    a = 1 ./ a; a[1] = 0

    ξ == 0 && return a
    α = fdiff_interpolation_expansion_coeffs(-ξ, k, fwd)

    a,α = Base.promote(a,α)

return polynom_product_expansion(a, α, k)

end
#...............................................................................
function bwd_differentiation_expansion_coeffs(ξ::T, k=3) where T<:Real
# ==============================================================================
#   Backward difference expansion coefficients for differentiation of an
#   analytic function f(x) tabulated under the convention f[n-k,n] and
#   evaluated at the interpolation position n+ξ.
# ==============================================================================
    Float = (Float64,Float32,Float16)

    b = T ∈ Float ? Base.ones(T,k+1) :
                    T <: Rational{}  ? Base.ones(T,k+1) :
                                       Base.ones(Rational{T},k+1)
    for i=1:k
        b[i+1] = i
    end

    b = 1 ./ b; b[1] = 0

    ξ == 0 && return b
    β = fdiff_interpolation_expansion_coeffs(ξ, k, bwd)

    b,β = Base.promote(b,β)

    return polynom_product_expansion(b, β, k)

end
#...............................................................................

@doc raw"""
    fdiff_differentiation_expansion_coeffs(ξ::T [, k=3 [, notation=bwd]]) where T<:Real

Finite-difference expansion coefficient vector defining ``k^{th}``-order
lagrangian *differentiation* of the tabulated analytic function ``f[n]``
at offset ``ξ`` (with respect to index position ``n``), which is positive
for increasing index and negative for decreasing index.

**Forward difference notation** (`notation = fwd`)

```math
\frac{df}{dξ}[n+ξ]=\sum_{p=0}^kα_p(ξ)Δ^{p}f[n]
```

Offset convention: ``ξ = -σ`` with respect to index ``n`` in tabulated
interval ``f[n:n+k]``

**Backward difference notation** (`notation = bwd`)

```math
\frac{df}{dξ}[n+ξ]=\sum_{p=0}^kβ_p(ξ)∇^{p}f[n]
```
where ``β(ξ) ≡ [β_0(ξ),\ ⋯,\ β_p(ξ)]``

Offset convention: ``ξ = σ`` with respect to index ``n`` in tabulated
interval ``f[n-k:n]``
#### Example:
```
k = 2; ξ = 0
o = fdiff_differentiation_expansion_coeffs(ξ, k); println(o)
 [0.0, 1.0, -1.5]
```
"""
function fdiff_differentiation_expansion_coeffs(ξ::T, k=3, notation=bwd) where T<:Real

    o = isforward(notation) ? fwd_differentiation_expansion_coeffs(ξ, k) :
                              bwd_differentiation_expansion_coeffs(ξ, k)
    return o

end

# ================= fdiff_differentiation(f, v; k=3) ==============================

@doc raw"""
    fdiff_differentiation(f::Vector{T}, v::V; k=3) where {T<:Real, V<:Real}

``k^{th}``-order (default *third* order) lagrangian *differentiation* of
the analytic function ``f``, tabulated in forward order on a uniform grid.
#### Example:
```
f = [x^3 for x=-5:5]; println(round.(Int,f))
l = length(f)
f′ = [fdiff_differentiation(f, v) for v=1:l]; println(round.(Int,f′))
  [-125, -64, -27, -8, -1, 0, 1, 8, 27, 64, 125]
  [75, 48, 27, 12, 3, 0, 3, 12, 27, 48, 75]

f′= fdiff_differentiation(f, 1) ; println("f′(1) = $(f′))
  f′(1) = 75//1

f′= fdiff_differentiation(f, 6.5) ; println("f′(6.5) = $(f′)")
    f′(6.5) = 0.75
```
For a cubic function the third-order lagrangian differentiation is exact -
see Figure below.

![Image](./assets/lagrangian_differentiation.png)
"""
function fdiff_differentiation(f::Vector{T}, v::V; k=3) where {T<:Real, V<:Real}

    l = length(f)
    k = min(k,l-1)
    k > 0 || error("Error: k ≥ 1 required for lagrangian interpolation")
    n = v ≤ 1 ? 1 : v < l-k ? floor(Int,v) : l-k
    ξ = v-n
    α = fdiff_differentiation_expansion_coeffs(ξ, k, fwd)
    w = fdiff_expansion_weights(α, fwd, reg)
    o = w ⋅ f[n:n+k]

    return o

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

    for i=0:k
        coeffs = CamiXon.fdiff_differentiation_expansion_coeffs(k-i,k)
        m[1+i,1:k+1] = fdiff_expansion_weights(coeffs)
    end

    return m

end

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
        a = CamiXon.fdiff_adams_moulton_expansion_coeffs(k)
        D = Base.denominator(Base.gcd(a))          # == Adams-Moulton devisor
        o = devisor ? (k, D, Base.round.(Int, o* D)) :
                      Base.round.(Int, o* D) // D  # convert to Rational{Int}
    end

    return o

end

# ======================== trapezoidal_integration(f, domain; k=5) =============

@doc raw"""
    trapezoidal_integration(f, x1, x2, weights)

Integral of the tabulated function ``f=[f_0,⋯\ f_n]`` over the `domain`
``x1 ≤ x ≤ x2`` using the optimized trapezoidal rule with endpoint correction
by the weights vector `weights`,
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
