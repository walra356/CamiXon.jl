# SPDX-License-Identifier: MIT

# author: Jook Walraven - 14-2-2023

# ==============================================================================
#                               grid_autoset.jl
# ==============================================================================# ....................... autoRmax(atom, orbit) ................................

@doc raw"""
    autoRmax(atom::Atom, orbit::Orbit)

Largest relevant radial distance in a.u. (rule of thumb value)
```math
    R_{max} = (2n^2 + 20n + 62)/Zc,
```
where ``n`` is the principal quantum number and ``Z_c`` the Rydberg charge
#### Example:
```
codata = castCodata(2018)
atom = castAtom(Z=1, A=1, Q=0)
orbit = castOrbit(n=1, ℓ=0)
rmax = autoRmax(atom::Atom, orbit::Orbit); println("rmax = $(rmax) a.u.")
    Element created: H, hydrogen, Z=1, weight=1.008
    Isotope created: ¹H, hydrogen, Z=1, A=1, N=0, R=0.8783, M=1.007825032, I=1/2⁺, μI=2.792847351, Q=0.0, RA=99.9855%, (stable)
    Atom created: hydrogen, neutral atom, ¹H, Z=1, A=1, Q=0, Zc=1
    Orbital: 1s
        principal quantum number: n = 1
        radial quantum number: n′ = 0 (number of nodes in radial wavefunction)
        orbital angular momentum of valence electron: ℓ = 0
    rmax = 63.0 a.u.
```
"""
function autoRmax(atom::Atom, orbit::Orbit)
# ==============================================================================
#  Discretization range in atomic units
# ==============================================================================
     n = orbit.n
     ℓ = orbit.ℓ
    Zc = atom.Zc

    #Rmax = 4(n^2+20)/Zc
    Rmax = (2n^2 + 20n + 62)/Zc
    #Rmax = (3n^2 -ℓ*(ℓ+1))/Zc

    return Rmax

end

# .......................... autoNtot(orbit) ...................................

@doc raw"""
    autoNtot(orbit::Orbit, Nboost=1)

Total number of gridpoints (rule of thumb value)
```math
    N_{tot} = (70 + 50 * n) * N_{boost},
```
where ``n`` is the principal quantum number and `Nboost` a multiplier to boost numerical precision
### Example:
```
orbit = castOrbit(n=1, ℓ=0)
autoNtot(orbit)
 Orbit created: 1s - (n = 1, n′ = 0, ℓ = 0)

 100
```
"""
function autoNtot(orbit::Orbit, Nboost=1)
# ==============================================================================
#  Total number of grid points
# ==============================================================================

    n = orbit.n

    Ntot = (70 + 50 * n) * Nboost

    return Ntot

end

# .................... autoPrecision(Rmax, orbit) ..............................

@doc raw"""
    autoPrecision(Rmax::T, orbit::Orbit) where T<:Real

Floating point precision (rule of thumb value)
### Example:
```
atom = castAtom(Z=1)
orbit = castOrbit(n=1,ℓ=0)
Rmax = autoRmax(atom, orbit)
o = autoPrecision(Rmax, orbit); println("precision = $o")
    Element created: H, hydrogen, Z=1, weight=1.008
    Isotope created: ¹H, hydrogen, Z=1, A=1, N=0, R=0.8783, M=1.007825032, I=1/2⁺, μI=2.792847351, Q=0.0, RA=99.9855%, (stable)
    Atom created: hydrogen, neutral atom, ¹H, Z=1, A=1, Q=0, Zc=1
    Orbital: 1s
        principal quantum number: n = 1
        radial quantum number: n′ = 0 (number of nodes in radial wavefunction)
        orbital angular momentum of valence electron: ℓ = 0
    precision = Float64
```
"""
function autoPrecision(Rmax::T, orbit::Orbit) where T<:Real
# ==============================================================================
#  Floating point precision (rule of thumb value)
# ==============================================================================

    ℓ = orbit.ℓ

    mytype = Rmax^(ℓ+1) == Inf ? BigFloat : T

    mytype ≠ T && println("\nautoPrecision: Rmax^(ℓ+1) => Inf - overflow protection: type promoted to BigFloat\n")

    return mytype

end

# ...................... autoSteps(Ntot, Rmax) .................................

@doc raw"""
    autoSteps(ID::Int, Ntot::Int, Rmax::T; p=5) where T<:Real
    autoSteps(ID::Int, Ntot::Int, Rmax::T; coords=[0,1]) where T<:Real

Step size parameter (h) and range parameter (r0) (rule of thumb values).
### Example:
```
(h, r0) = autoSteps(1, 100, 100)
    (0.1, 0.004540199100968777)
```
"""
function autoSteps(ID::Int, Ntot::Int, Rmax::T; p=5, coords=[0,1]) where T<:Real
# ==============================================================================
#  Step size parameter (h) and range parameter (r0)
# ==============================================================================

    h = 10/Ntot

    r0 = Rmax / gridfunction(ID, Ntot, h; p, coords)

    return h, r0

end

# ..............................................................................

@doc raw"""
    gridfunction(ID::Int, n::Int, h::T; p=5, coords=[0,1], deriv=0) where T <: Real

* `ID = 1`: exponential grid function,
```math
    f[n] = \text{exp}(h(n-1)) - 1.0
```
* `ID = 2`: quasi-exponential grid function degree `p` (linear grid for `p = 1`),
```math
    f[n] = h(n-1) + \frac{1}{2}(h(n-1))^2 + ⋯ + \frac{1}{p!}(h(n-1))^p
```
* `ID = 3`: linear grid function,
```math
    f[n] = h(n-1)
```
* `ID = 4`: polynomial grid function of degree `p = length(c)` based on `polynom` ``c = [c_1,c_2,⋯\ c_p]``,
```math
    f[n] = c_1h(n-1) + c_2(h(n-1))^2 + ⋯ + c_p(h(n-1))^p
```
#### Examples:
```
h = 0.1
r = [gridfunction(1, n-1, h) for n=1:5]                            # exponential
 [0.0, 0.10517091807564771, 0.22140275816016985, 0.3498588075760032, 0.49182469764127035]

r = [gridfunction(2, n-1, h; p = 4) for n=1:5]  # quasi exponential (degree p=4)
 [0.0, 0.10517083333333321, 0.22140000000000004, 0.3498375, 0.49173333333333336]

r = [gridfunction(3, n-1, h) for n=1:5]              # linear
  [0.0, 0.1, 0.2, 0.3, 0.4]

r′= [gridfunction(3, n-1, h; deriv=1) for n=1:5]     # linear (first derivative)
   [0.1, 0.1, 0.1, 0.1, 0.1]

  r = [gridfunction(4, n-1, h; coords = [0,1,1/2,1/6,1/24]) for n=1:5]  # polynomial of degree 4)
   [0.0, 0.10517083333333334, 0.2214, 0.3498375000000001, 0.49173333333333336]
```
"""
function gridfunction(ID::Int, n::Int, h::T; p=5, coords=[0,1], deriv=0) where T <: Real

    ID == 1 && return _walterjohnson(n, h; deriv)
    ID == 2 && return _jw_gridfunction(n, h; deriv, p)
    ID == 3 && return _linear_gridfunction(n, h; deriv)
    ID == 4 && return CamiMath.polynomial(coords, h*n; deriv)

    return error("Error: unknown gridfunction")

end

# = autoGrid(atom, orbit, T; p=0, coords=[], Nboost=1, epn=5, k=7, msg=true) ===

@doc raw"""
    autoGrid(atom, orbit,  T; Nboost=1, epn=5, k=7, msg=true, p=0)
    autoGrid(atom, orbits, T; Nboost=1, epn=5, k=7, msg=true, p=0)
    autoGrid(atom, orbit,  T; Nboost=1, epn=5, k=7, msg=true, coords=[])
    autoGrid(atom, orbits, T; Nboost=1, epn=5, k=7, msg=true, coords=[])

Automatic setting of grid parameters for a given orbit [`Orbit`](@ref) or an
array of orbits - `orbits = [orbit1, orbit2, ⋯]`. Important cases:
* `p == 0` (exponential radial grid)
* `p == 1` (linear radial grid)
* `p > 1` (quasi-exponential radial grid)
* `coords=[]` (free polynomial grid based on the `coords`)
* `Nboost` (multiplier to boost numerical precision)
* `epn` (endpoint number: odd number to be used for trapezoidal integration with endpoint correction)
* `k` (Adams-Moulton order to be used for `k+1`-point Adams-Moulton integration)
#### Example:
```
codata = castCodata(2018)
atom = castAtom(;Z=1, A=1, Q=0, msg=false)
orbit = castOrbit(n=75, ℓ=0, msg=false)
grid = autoGrid(atom, orbit, Float64);
    Grid created: exponential, Float64, Rmax = 16935.0 a.u., Ntot = 3800, h = 0.00263158, r0 = 0.768883

plot_gridfunction(grid, 1:grid.N; title="")
```
The plot is made using CairomMakie.
NB.: `plot_gridfunction` is not part of the `CamiXon` package.
![Image](./assets/exponential_grid.png)
"""
function autoGrid(atom::Atom, orbit::Orbit, T::Type; p=0, coords=[], Nboost=1, epn=5, k=7, msg=false)

    T ∈ [Float64,BigFloat] || println("autoGrid: grid.T = $T => Float64 (was enforced by automatic type promotion)")

    ID = (p < 1) & (length(coords) < 2) ? 1 :
         (p ≥ 1) & (length(coords) < 2) ? (p == 1 ? 3 : 2) :
         (p < 1) & (length(coords) ≥ 2) ? 4 : error("Error: unknown grid")

    Ntot = autoNtot(orbit, Nboost)
    Rmax = autoRmax(atom, orbit)

    T = T == BigFloat ? T : autoPrecision(Rmax, orbit)
    h, r0 = autoSteps(ID, Ntot, Rmax; p, coords)

    return castGrid(ID, Ntot, T; h, r0, p, coords, epn, k, msg)

end
function autoGrid(atom::Atom, orbits::Vector{Orbit}, T::Type; p=0, coords=[], Nboost=1, epn=5, k=7, msg=false)

    T ∈ [Float64,BigFloat] || println("autoGrid: grid.T = $T => Float64 (enforced by automatic type promotion)")

    ID = (p < 1) & (length(coords) < 2) ? 1 :
         (p ≥ 1) & (length(coords) < 2) ? (p == 1 ? 4 : 2) :
         (p < 1) & (length(coords) ≥ 2) ? 3 : error("Error: unknown grid")

    R = round.(Int, [autoRmax(atom, orbits[i]) for i ∈ eachindex(orbits)])
    Rmax = maximum(R)
    Ntot = maximum([autoNtot(orbits[i], Nboost)*Rmax÷R[i] for i ∈ eachindex(orbits)])
    h,r0 = autoSteps(ID, Ntot, Rmax; p, coords)

    for i ∈ eachindex(orbits)
        autoPrecision(Rmax, orbits[i]) == BigFloat ? T = BigFloat : false
    end

    return castGrid(ID, Ntot, T; h, r0, p, coords, epn, k, msg)

end
