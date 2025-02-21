# SPDX-License-Identifier: MIT

# author: Jook Walraven - 14-2-2023

# ==============================================================================
#                               grid_autoset.jl
# ==============================================================================

# ....................... autoRmax(atom, orbit) ................................

@doc raw"""
    autoRmax(rmax::T, atom::Atom, orbit::Orbit) where T<:Real

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
rmax = autoRmax(rmax, atom, orbit); println("rmax = $(rmax) a.u.")
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
function autoRmax(atom::Atom, orbit::Orbit; rmax=0.0)
# ==============================================================================
#   Discretization range in atomic units
# ==============================================================================
    n = orbit.n
    ℓ = orbit.ℓ
    Zc = atom.Zc
        
    rmax = rmax > 0 ? rmax : (2.5n^2 + 75.0)/Zc
    
    return rmax
        
end

# .......................... autoNtot(orbit) ...................................

@doc raw"""
    autoNtot(orbit::Orbit)

Total number of gridpoints (rule of thumb value)
```math
    N_{tot} = 400 + 100 n,
```
where ``n`` is the principal quantum number
### Example:
```
orbit = castOrbit(n=1, ℓ=0)
autoNtot(orbit)
 Orbit created: 1s - (n = 1, n′ = 0, ℓ = 0)

 500
```
"""
function autoNtot(orbit::Orbit)
# ==============================================================================
#   Total number of grid points
# ==============================================================================
    N = 400 + 100 * (orbit.n)
        
    return N
    
end
    
# .................... autoPrecision(rmax, orbit) ..............................

@doc raw"""
    autoPrecision(rmax::T, orbit::Orbit) where T<:Real 

Floating point precision (rule of thumb value)
### Example:
```
atom = castAtom(Z=1)
orbit = castOrbit(n=1,ℓ=0)
rmax = autoRmax(rmax, atom, orbit)
o = autoPrecision(rmax, orbit); println("precision = $o")
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
function autoPrecision(rmax::T, orbit::Orbit) where T<:Real
# ==============================================================================
#  Floating point precision (rule of thumb value)
# ==============================================================================

    ℓ = orbit.ℓ

    mytype = rmax^(ℓ+1) == Inf ? BigFloat : Float64

    mytype ≠ T && println("\nautoPrecision: rmax^(ℓ+1) => Inf - overflow protection: type promoted to BigFloat\n")

    return mytype

end

@doc raw"""
    autoGrid(atom, orbit,  T; Nboost=1, epn=5, k=7, msg=true, p=0)
    autoGrid(atom, orbits, T; Nboost=1, epn=5, k=7, msg=true, p=0)
    autoGrid(atom, orbit,  T; Nboost=1, epn=5, k=7, msg=true, polynom=[])
    autoGrid(atom, orbits, T; Nboost=1, epn=5, k=7, msg=true, polynom=[])

Automatic setting of grid parameters for a given orbit [`Orbit`](@ref) or an
array of orbits - `orbits = [orbit1, orbit2, ⋯]`. Important cases:
* `p == 0` (exponential radial grid)
* `p == 1` (linear radial grid)
* `p > 1` (quasi-exponential radial grid)
* `polynom=[]` (free polynomial grid based on the `polynom`)
* `Nboost` (multiplier to boost numerical precision)
* `epn` (endpoint number: odd number to be used for trapezoidal integration with endpoint correction)
* `k` (Adams-Moulton order to be used for `k+1`-point Adams-Moulton integration)
#### Example:
```
codata = castCodata(2018)
atom = castAtom(;Z=1, A=1, Q=0, msg=false)
orbit = castOrbit(n=75, ℓ=0, msg=false)
grid = autoGrid(atom, orbit, Float64);
    CamiDiff.Grid created: exponential, Float64, rmax = 16935.0 a.u., N = 3800, h = 0.00263158, r0 = 0.768883

plot_gridfunction(grid, 1:grid.N; title="")
```
The plot is made using CairomMakie.
NB.: `plot_gridfunction` is not part of the `CamiXon` package.
![Image](./assets/exponential_grid.png)
"""
function autoGrid(atom::Atom, orbit::Orbit, T::Type; h=0, p=0, polynom=[], N=0, rmax=0, epn=5, k=5, msg=false)

    T ∈ [Float64,BigFloat] || println("autoGrid: grid.T = $T => Float64 (was enforced by automatic type promotion)")

    ID = (p < 1) & (length(polynom) < 2) ? 1 :
         (p ≥ 1) & (length(polynom) < 2) ? (p == 1 ? 3 : 2) :
         (p < 1) & (length(polynom) ≥ 2) ? 4 : error("Error: unknown grid")

    N = N == 0 ? autoNtot(orbit) : N
    rmax = autoRmax(atom, orbit; rmax)

    #T = T == BigFloat ? T : autoPrecision(rmax, orbit)
    h = h ≠ 0 ? h : N < 100 ? T(1//N) : T(1//100)
    h = T(10//N)

    return CamiDiff.castGrid(ID, N, T; h, rmax, p, polynom, epn, k, msg)

end