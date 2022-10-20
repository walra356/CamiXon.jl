# ============================ hydrogen sector =================================

# ..............................................................................
function _hydrogenic_norm(n::U, ℓ::U) where U<:Real

    T = (n + ℓ) > 20 ? BigInt : Int

    o = T(1)

    for i=(n-ℓ):(n+ℓ)
        o *= T(i)
    end

    o *= T(2n)

    return o

end
# ..............................................................................

@doc raw"""
    hydrogenic_reduced_wavefunction(atom::Atom, orbit::Orbit, grid::Grid, def::Def)

Analytic expression for the hydrogenic wavefunction written in the format
```math
    Z = \tilde{χ} + im \tilde{χ}^′``,
```
where ``\tilde{χ}_{nℓ}(ρ)`` is the *reduced* radial wavefunction and
``\tilde{χ}^′_{nℓ}(ρ)`` its derivative, with ``ρ`` the radial distance
to the nucleus in a.u.. The expression is evaluated for a given [`Atom`](@ref)
in a given [`Orbit`](@ref) on a given [`Grid`](@ref). The argument [`Def`](@ref)
completes the definition of the problem.
```math
    \tilde{\chi}_{nl}(\rho)
    =\mathcal{N}_{nl}^{-1/2}(2Z/n)^{l+3/2}\rho^{l+1}e^{-Zρ/n}
    L_{n-l-1}^{2l+1}(2Z\rho/n)
```
where ``L_{n-l-1}^{2l+1}(2Z\rho/n)`` is the generalized Laguerre polynomial
[`generalized_laguerreL`](@ref) and
```math
    \mathcal{N}_{nl}
    = {\displaystyle \int\nolimits _{0}^{\infty}}x^{2l+2}e^{-x}
    \left[L_{n-l-1}^{2l+1}(x)\right]^{2}dx
    = \frac{2n\Gamma(n+l+1)}{\Gamma(n-l)}
```
is the norm of the wavefunction.
#### Example:
```
atom = castAtom(Z=1, A=1, Q=0)
orbit = castOrbit(n=25, ℓ=10)
grid = autoGrid(atom, orbit, Float64; Nboost=1, msg=true)
def = castDef(grid, atom, orbit, codata)
Z = hydrogenic_reduced_wavefunction(atom, orbit, grid, def);
    Element created: H, hydrogen, Z=1, weight=1.008
    Isotope created: ¹H, hydrogen, Z=1, A=1, N=0, R=0.8783, M=1.007825032, I=1/2⁺, μI=2.792847351, Q=0.0, RA=99.9855%, (stable)
    Atom created: hydrogen, neutral atom, ¹H, Z=1, A=1, Q=0, Zc=1
    Orbital: 25n
    principal quantum number: n = 25
    radial quantum number: n′ = 14 (number of nodes in radial wavefunction)
    orbital angular momentum of valence electron: ℓ = 10
    Grid created: exponential, Float64, Rmax = 1935.0 a.u., Ntot = 1300, h = 0.00769231, r0 = 0.0878529
    Def created for hydrogen 25n on exponential grid of 1300 points

plot_wavefunction(Z, 1:grid.N, grid, def)
```
The plot is made using `CairomMakie`.
NB.: `plot_wavefunction` is not included in the `CamiXon` package.
![Image](./assets/H1_25n.png)
"""
function hydrogenic_reduced_wavefunction(atom::Atom, orbit::Orbit, grid::Grid, def::Def)

    Zval = atom.Z
    n = orbit.n
    ℓ = orbit.ℓ
    r = grid.r

    (atom.Z - atom.Q) == 1 || error("Error: Z-Q ≠ 1 (atom not hydrogenic)")

    grid.N == def.pos.N || error("Error: grid.N ≠ def.pos.N")

    norm = _hydrogenic_norm(n, ℓ)

    a = float(big(2Zval)//big(n))
    b = a^(ℓ+1)*sqrt(a/big(norm))

    coords = float(generalized_laguerre_coords(n-ℓ-1, 2ℓ+1))

    P = b .* [r[i]^(ℓ+1) * exp(-0.5a*r[i]) * polynomial(coords, a*r[i]) for i ∈ eachindex(r)]
    Q = b .* [r[i]^ℓ * exp(-0.5a*r[i]) * (((ℓ+1)-0.5a*r[i]) * polynomial(coords, a*r[i]) + a*r[i]*polynomial(coords, a*r[i]; deriv=1)) for i ∈ eachindex(r)]

    return P + im * Q

end

# ======================== bohrformula(atom, term) =============================

@doc raw"""
    bohrformula(Z::Int, n::Int)

Hydrogenic energy (in Hartree a.u.) for *atom* with *atomic number* `Z` and
*principal quantum number* `n`.
```math
    E_n = - \frac{Z^2}{2n^2}
```
#### Example:
```
Z = 2
n = 4
bohrformula(Z,n)
 -0.125
```
"""
bohrformula(Z::Int, n::Int) = -(1//2)*(Z//n)^2

# ========================= demo_hydrogen(; n=3, ℓ=2) ==========================

@doc raw"""
    demo_hydrogen(; n=3, ℓ=2, codata=castCodata(2018))

Solves Schrödinger equation for hydrogen atom with principal quantum number `n`
and rotational quantum number `ℓ`.

#### Example:
```
Ecal, grid, def, adams = demo_hydrogen(n=1, ℓ=0);
    Def created for hydrogen 1s on exponential grid of 100 points

E = 1.5Ecal
E, def, adams, Z = adams_moulton_master(E, grid, def, adams; Δν=Value(1,"kHz"), imax=25, msg=true);

plot_wavefunction(Z, 1:def.pos.N, grid, def; reduced=false)
```
The plot is made using `CairomMakie`. Note the discontinuity in the derivative.
NB.: `plot_wavefunction` is not included in the `CamiXon` package.
![Image](./assets/hydrogen-1s.png)
"""
function demo_hydrogen(; n=3, ℓ=2, codata=castCodata(2018))

    atom = castAtom(;Z=1, A=1, Q=0, msg=false)
    orbit = castOrbit(; n, ℓ, msg=false)
    grid = autoGrid(atom, orbit, Float64; msg=false)
    def = castDef(grid, atom, orbit, codata; msg=true )
    E = convert(grid.T, bohrformula(atom.Z, orbit.n))
    adams = castAdams(E, grid, def)

    return E, grid, def, adams

end

# =================== restore_wavefunction(Z, grid) ============================

@doc raw"""
    restore_wavefunction(Z::Vector{Complex{T}}, grid::Grid{V}) where {T<:Real, V<:Real}

Conversion from the *reduced* radial wavefunction ``\tilde{\chi}_{nl}(ρ)``
to the *ordinary* radial wavefuntion ``\tilde{R}_{nl}(ρ)``,
```math
    \tilde{R}_{nl}(ρ)=\tilde{\chi}_{nl}(ρ)/ρ,
```
where ``ρ`` is the radial distance to the nucleus in a.u..
#### Example:
```
atom = castAtom(Z=1, A=1, Q=0; msg=false);
orbit = castOrbit(n=1, ℓ=0; msg=false);
grid = autoGrid(atom, orbit, Float64; Nboost=1, msg=false);
def = castDef(grid, atom, orbit, codata);

RH1s_example = [RH1s(atom.Z, grid.r[n]) for n=1:grid.N]
XH1s_generic = hydrogenic_reduced_wavefunction(atom, orbit, grid, def)

XH1s_example = reduce_wavefunction(RH1s_example)
RH1s_generic = restore_wavefunction(XH1s_generic)

XH1s_example ≈ XH1s_generic
    true

RH1s_example ≈ RH1s_generic
    true
```
"""
function restore_wavefunction(Z::Vector{Complex{T}}, grid::Grid{V}) where {T<:Real, V<:Real}

    χ = real(Z)
    χ′= imag(Z)
    r = grid.r

    R = χ ./ r
    R′= (χ′ .- χ ./ r) ./ r

    R[1] = fdiff_interpolation(R[2:end], 0) # extrapolate to r=0 to handle division by "zero"
    R′[1] = fdiff_interpolation(R′[2:end], 0)


    return R + im * R′

end
# =================== reduce_wavefunction(Z, grid) ============================

@doc raw"""
    reduce_wavefunction(Z::Vector{Complex{T}}, grid::Grid{V}) where {T<:Real, V<:Real}

Conversion from the *ordinary* radial wavefunction ``\tilde{R}_{nl}(ρ)``
to the *reduced* radial wavefuntion
```math
    \tilde{\chi}_{nl}(ρ) = ρ \tilde{R}_{nl}(ρ).
```
where ``ρ`` is the radial distance to the nucleus in a.u..
#### Example:
```
atom = castAtom(Z=1, A=1, Q=0; msg=false);
orbit = castOrbit(n=1, ℓ=0; msg=false);
grid = autoGrid(atom, orbit, Float64; Nboost=1, msg=false);
def = castDef(grid, atom, orbit, codata);

RH1s_example = [RH1s(atom.Z, grid.r[n]) for n=1:grid.N]
XH1s_generic = hydrogenic_reduced_wavefunction(atom, orbit, grid, def)

XH1s_example = reduce_wavefunction(RH1s_example)
RH1s_generic = restore_wavefunction(XH1s_generic)

XH1s_example ≈ XH1s_generic
    true

RH1s_example ≈ RH1s_generic
    true
```
"""
function reduce_wavefunction(Z::Vector{Complex{T}}, grid::Grid{V}) where {T<:Real, V<:Real}

    R = real(Z)
    R′= imag(Z)
    r = grid.r

    χ = r .* R
    χ′= r .* R′ + R

    return χ + im * χ′

end

# =============================== RHqs(Z, r) ===================================

@doc raw"""
    RH1s(Z::Int, r::T) where T <:Real

Analytic expression for the *hydrogenic* 1s radial wavefunction
and its derivative in the format,
```math
    Z = \tilde{R} + i \tilde{R}^′,
```
where
```math
    \tilde{R}_{1s}(ρ) = Z^{3/2} 2 e^{-Zρ}
```
is the radial wavefunction and
```math
    \tilde{R}^′_{1s}(ρ) = -Z^{5/2} 2 e^{-Zρ}
```
its derivative, with ``ρ`` the radial distance to the nucleus in a.u..
#### Example:
```
atom = castAtom(Z=1, A=1, Q=0; msg=false);
orbit = castOrbit(n=1, ℓ=0; msg=false);
grid = autoGrid(atom, orbit, Float64; Nboost=1, msg=false);
def = castDef(grid, atom, orbit, codata);

RH1s_example = [RH1s(atom.Z, grid.r[n]) for n=1:grid.N];

plot_wavefunction(RH1s_example, 1:grid.N, grid, def; reduced=false)
```
The plot is made using `CairomMakie`.
NB.: `plot_function` is not included in the `CamiXon` package.
![Image](./assets/RH1s.png)
"""
function RH1s(Z::Int, r::T) where T <:Real

    p = 2.0 * Z^(3/2) * exp(-Z*r)
    q = -2.0 * Z^(5/2) * exp(-Z*r)

    return p + im * q

end

# =============================== RH2p(Z, r) ===================================

@doc raw"""
    RH2p(Z::Int, r::T) where T <:Real

Analytic expression for the *hydrogenic* 1s reduced radial wavefunction
and its derivative in a complex number format,
```math
    \tilde{R}_{2p}(\rho)&=\left(Z/2\right)^{3/2}\sqrt{1/3}(Z\rho/2)2ρe^{-Zρ/2}
```
#### Example:
```
atom = castAtom(Z=1, A=1, Q=0; msg=false);
orbit = castOrbit(n=2, ℓ=1; msg=false);
grid = autoGrid(atom, orbit, Float64; Nboost=1, msg=false);
def = castDef(grid, atom, orbit, codata)

RH2p_example = [RH2p(atom.Z, grid.r[n]) for n=1:grid.N]

plot_wavefunction(RH2p_example, 1:grid.N, grid, def; reduced=false)
```
The plot is made using `CairomMakie`.
NB.: `plot_wavefunction` is not included in the `CamiXon` package.
![Image](./assets/RH2p.png)
"""
function RH2p(Z::Int, r::T) where T <:Real

    P = 2.0 * (Z/2)^(3/2) * sqrt(1/3) * (Z*r/2.0) * exp(-Z*r/2.0)
    Q = 2.0 * (Z/2)^(5/2) * sqrt(1/3) * exp(-Z*r/2.0) * (1.0 - Z*r/2.0)

    return P + im * Q

end

χHe1s(r) = 4.0 * sqrt(2) * r * exp(-2.0r) + im * 4.0 * sqrt(2) * exp(-2.0r) * (1 - 2.0r)
gridHe1s(grid) = [χHe1s(grid.r[n]) for n=1:grid.N]
