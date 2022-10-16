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
    hydrogenic_wavefunction(atom::Atom, orbit::Orbit, grid::Grid, def::Def)

Analytic function for the hydrogenic *reduced* radial wavefunction of a given
[`atom`](@ref) in a given [`orbit`](@ref) on a given [`grid`](@ref).
The argument [`def`](@ref) completes the definition of the problem.

```math
    \tilde{\chi}_{nl}(\rho)
    =\mathcal{N}_{nl}^{-1/2}(2Z/n)^{l+3/2}\rho^{l+1}e^{-Z\rho/n}
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
Z = hydrogenic_wavefunction(atom, orbit, grid, def);
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
function hydrogenic_wavefunction(atom::Atom, orbit::Orbit, grid::Grid, def::Def)

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
The plot is made using CairomMakie. Note the discontinuity in the derivative.
NB.: `plot_wavefunction` is not included in the `CamiXon` package.
```
Ecal, grid, def, adams = demo_hydrogen(n=1, ℓ=0);
    Def created for hydrogen 1s on exponential grid of 100 points

E = 1.5Ecal
E, def, adams, Z = adams_moulton_master(E, grid, def, adams; Δν=Value(1,"kHz"), imax=25, msg=true);

plot_wavefunction(Z, 1:def.pos.N, grid, def; undo_reduction=true)
```
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
