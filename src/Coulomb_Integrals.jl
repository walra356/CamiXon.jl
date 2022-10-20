@doc raw"""
    a_coeff(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)

Angular coefficient for the *direct* Coulomb integral:

```math
a^{k}(lm_{l};l^{\prime}m_{l^{\prime}})=(-)^{m_{l}+m_{l^{\prime}}}
(2l+1)(2l^{\prime}+1)\left(\begin{array}{ccc}
l & k & l\\
0 & 0 & 0
\end{array}\right)\left(\begin{array}{ccc}
l & k & l\\
-m_{l} & 0 & m_{l}
\end{array}\right)\left(\begin{array}{ccc}
l^{\prime} & k & l^{\prime}\\
0 & 0 & 0
\end{array}\right)\left(\begin{array}{ccc}
l^{\prime} & k & l^{\prime}\\
-m_{l^{\prime}} & 0 & m_{l^{\prime}}
\end{array}\right)
```
#### Example:
```
a(2,1,1,2,2)
    2//35

a(6,3,2,3,-1)
    -250//20449
```
"""
function a_coeff(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)

    Base.iseven(k) || return 0
    0 ≤ k ≤ 2min(l,l′) || return 0
    k > 0 || return 1

    Δ1 = triangle_coefficient(l, k, l)
    Δ2 = triangle_coefficient(l′, k, l′)
    T = _threeJroot2(l, 0, k, 0, l, 0) * _threeJroot2(l, -ml, k, 0, l, ml) * _threeJroot2(l′, 0, k, 0, l′, 0) * _threeJroot2(l′, -ml′, k, 0, l′, ml′)
    S = _Racah_sum(l, 0, k, 0, l) * _Racah_sum(l, -ml, k, 0, l) * _Racah_sum(l′, 0, k, 0, l′) * _Racah_sum(l′, -ml′, k, 0, l′)

    a = (2l+1) * (2l′+1) * Δ1 * Δ2 * S
    o = a * a * T
    num = Int(sqrt(numerator(o)))
    den = Int(sqrt(denominator(o)))

    o = sign(S) * num//den

    return o

end

@doc raw"""
    b_coeff(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)

Angular coefficient for the *exchange* Coulomb integral:

```math
b^{k}(lm_{l};l^{\prime}m_{l^{\prime}})=(2l+1)(2l^{\prime}+1)
\left(\begin{array}{ccc}
l & k & l^{\prime}\\
0 & 0 & 0
\end{array}\right)^{2}\left(\begin{array}{ccc}
l & k & l^{\prime}\\
-m_{l} & (m_{l}-m_{l^{\prime}}) & m_{l^{\prime}}
\end{array}\right)^{2}
```
#### Example:
```
b(1,1,1,2,2)
    2//5

b(6,3,2,3,-1)
    1050//20449
```
"""
function b_coeff(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)

    Base.iseven(k+l+l′) || return 0
    abs(l-l′) ≤ k ≤ l+l′ || return 0

    Δ = triangle_coefficient(l, k, l′)
    T = _threeJroot2(l, 0, k, 0, l′, 0) * _threeJroot2(l, -ml, k, ml-ml′, l′, ml′)
    S = _Racah_sum(l, 0, k, 0, l′) * _Racah_sum(l, -ml, k, ml-ml′, l′)

    o = abs((2l+1) * (2l′+1) * Δ * Δ * S * S * T)
    num = Int(numerator(o))
    den = Int(denominator(o))

    o = num//den

    return o

end

# ======================== potUG(k, Z, grid) ===================================
@doc raw"""
    potUG(k::Int, Z1::Vector{Complex{T}}, Z1::Vector{Complex{T}}, grid::Grid{V}) where {T<:Real, V<:Real}

Coulomb integral for *exchange* screening,

```math
    U_{G}^{k}(\rho)
    =\frac{1}{\rho^{k+1}}\int_{0}^{\rho}\varrho^{k}\tilde{R}_{nl}(\varrho)
    \tilde{R}_{n^{\prime}l^{\prime}}(\varrho)\,
    \varrho^{2}d\varrho+\rho^{k}\int_{\rho}^{\infty}
    \frac{1}{\varrho^{k+1}}\tilde{R}_{nl}(\varrho)
    \tilde{R}_{n^{\prime}l^{\prime}}(\varrho)\,\varrho^{2}d\varrho.
```
#### Example:
```
atom = castAtom(Z=2, A=4, Q=0; msg=true)
orbit1 = castOrbit(n=1, ℓ=0; msg=true)
orbit2 = castOrbit(n=2, ℓ=0; msg=true)
scr = nothing
grid = autoGrid(atom, [orbit1,orbit2], Float64; Nboost=1, msg=true)
def1 = castDef(grid, atom, orbit1, codata; scr)
E = initE(def1)
adams = castAdams(E, grid, def1)
E, def, adams, Z1 = adams_moulton_master(E, grid, def1, adams; Δν=Value(1,"kHz"), imax=50, msg=false);

def2 = castDef(grid, atom, orbit2, codata; scr)
E = initE(def2)
adams = castAdams(E, grid, def2)
E, def, adams, Z2 = adams_moulton_master(E, grid, def2, adams; Δν=Value(1,"kHz"), imax=50, msg=false);

f = potUG(0, Z1, Z2, grid);
plot_function(f, 1:grid.N, grid; title="He4(1s,2s):  exchange screening potential")
```
The plot is made using `CairomMakie`.
NB.: `plot_function` is not included in the `CamiXon` package.
![Image](./assets/He41s-UG0.png)
"""
function potUG(k::Int, Z1::Vector{Complex{T}}, Z2::Vector{Complex{T}}, grid::Grid{V}) where {T<:Real, V<:Real}

    N = grid.N
    r = grid.r

    potUG_inner = [grid_integration(r.^k .* real(Z1) .* real(Z2), 1:n, grid) for n=2:N]
    potUG_outer = [grid_integration((1.0 ./ r).^(k+1) .* real(Z1) .* real(Z2), n:N, grid) for n=2:N]

    o = (potUG_inner .* r[2:N].^-(k+1)) .+ (potUG_outer .* r[2:N].^k)

    p = fdiff_interpolation(o, 0; k=4)

    prepend!(o,p)

    return o

end

# ======================== potUF(k, Z, grid) ===================================
@doc raw"""
    potUF(k::Int, Z::Vector{Complex{T}}, grid::Grid{V}) where {T<:Real, V<:Real}

Coulomb integral for *directe* screening,

```math
    U_{F}^{k}(\rho)
    =\frac{1}{\rho^{k+1}}\int_{0}^{\rho}\varrho^{k}
    \left[\tilde{R}_{nl}(\varrho)\right]^{2}
    \varrho^{2}d\varrho+\rho^{k}\int_{\rho}^{\infty}
    \frac{1}{\varrho^{k+1}}
    \left[\tilde{R}_{nl}(\varrho)\right]^{2}\varrho^{2}d\varrho.
```
#### Example:
```
atom = castAtom(Z=2, A=4, Q=0; msg=true)
orbit1 = castOrbit(n=1, ℓ=0; msg=true)
scr = nothing
grid = autoGrid(atom, orbit1, Float64; Nboost=1, msg=true)
def1 = castDef(grid, atom, orbit1, codata; scr)
E = initE(def1)
adams = castAdams(E, grid, def1)
E, def, adams, Z1 = adams_moulton_master(E, grid, def1, adams; Δν=Value(1,"kHz"), imax=50, msg=false);

pot = potUF(0, Z1, grid);
plot_function(pot, 1:grid.N, grid; title="He4(1s,1s):  direct screening potential")
```
The plot is made using `CairomMakie`.
NB.: `plot_function` is not included in the `CamiXon` package.
![Image](./assets/He41s-UF0.png)
"""
function potUF(k::Int, Z::Vector{Complex{T}}, grid::Grid{V}) where {T<:Real, V<:Real}

    return potUG(k, Z, Z, grid)

end
