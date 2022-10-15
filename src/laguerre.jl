# ============================= laguerre_coords(n::Int, α::T) where T<:Real =====================

@doc raw"""
    generalized_laguerre_coords(n::Int, α::T) where T<:Real

The coefficients of the generalized Laguerre polynomals of degree `n` for
parameter `α`.

```math
    c(n, α)[m] = \frac{\Gamma(α+n+1)}{\Gamma(α+m+1)}
    \frac{(-1)^{m}}{(n-m)!}\frac{1}{m!}
```
#### Example:
```
o = generalized_laguerre_coords(8,3); println(o)
    Rational{Int64}[165//1, -330//1, 231//1, -77//1, 55//4, -11//8, 11//144, -11//5040, 1//40320]
```
"""
function generalized_laguerre_coords(n::Int, α::T) where T<:Real

    coords = [_generalized_laguerre_coord(n, α, m) for m=0:n]

    return coords

end

# ============================= laguerre_coords(n::Int, α::T) where T<:Real =====================

@doc raw"""
    laguerre_coords(n::Int) where T<:Real

The coefficients of the Laguerre polynomals of degree `n`.

```math
    c(n)[m] = \frac{\Gamma(n+1)}{\Gamma(m+1)}\frac{(-1)^{m}}{(n-m)!}\frac{1}{m!}
```
#### Example:
```
o = laguerre_coords(8); println(o)
    Rational{Int64}[1//1, -8//1, 14//1, -28//3, 35//12, -7//15, 7//180, -1//630, 1//40320]
```
"""
function laguerre_coords(n::T) where T<:Real

    coords = [_generalized_laguerre_coord(n, 0, m) for m=0:n]

    return coords

end

@doc raw"""
    generalized_laguerreL(n::Int, α::U, x::V) where {U<:Real, V<:Real}

Generalized Laguerre polynomal of degree `n` for parameter `α`,

```math
    L_{n}^{α}(x)
    = \frac{1}{n!}e^{x}x^{-α}\frac{d^{n}}{dx^{n}}(e^{-x}x^{n+α})
    = \sum_{m=0}^{n}(-1)^{m}\binom{n+α}{n-m}\frac{x^{m}}{m!}
    = \sum_{m=0}^{n}c(n,α)[m]x^{m}
```
where ``c(n,α)[m]`` is the generalized Laguerre coordinate from
[`generalized_laguerre_coords`](@ref).
#### Example:
```
(xmin, Δx, xmax) = (0, 0.1, 11)
n = 8
α = -0.3
gL = [generalized_laguerreL(n, α, x) for x=xmin:Δx:xmax]
f = Float64.(gL);

plot_function(f, xmin, Δx, xmax; title="Laguere polynomial (of degree $n for α =$α)")
```
The plot is made using `CairomMakie`.
NB.: `plot_function` is not included in the `CamiXon` package.

![Image](./assets/laguerreL8.png)
"""
function generalized_laguerreL(n::Int, α::U, x::V) where {U<:Real, V<:Real}

    coords = convert.(BigFloat,generalized_laguerre_coords(n, α) )

    o = polynomial(coords, x; deriv=0)

    return o

end

@doc raw"""
    laguerreL(n::Int, x::T) where T<:Real

Laguerre polynomal of degree `n`,

```math
    L_{n}(x)
    = \frac{1}{n!}e^{x}\frac{d^{n}}{dx^{n}}(e^{-x}x^{n})
    = \sum_{m=0}^{n}(-1)^{m}\binom{n}{n-m}\frac{x^{m}}{m!}
    = \sum_{m=0}^{n}c(n)[m]x^{m}
```
where ``c(n)[m]`` is the Laguerre coordinate from [`laguerre_coords`](@ref).
#### Example:
```
(xmin, Δx, xmax) = (0, 0.1, 11)
n = 8
L = [laguerreL(n, x) for x=xmin:Δx:xmax]
f = Float64.(L);

plot_function(f, xmin, Δx, xmax; title="Laguere polynomial (of degree $n)")
```
The plot is made using `CairomMakie`.
NB.: `plot_function` is not included in the `CamiXon` package.
![Image](./assets/laguerreL8.png)
"""
function laguerreL(n::Int, x::T) where T<:Real

    coords = convert.(BigFloat,generalized_laguerre_coords(n, 0) )

    o = polynomial(coords, x; deriv=0)

    return o

end
