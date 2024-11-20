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
function autoRmax(atom::Atom, orbit::Orbit) # kanweg
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
function autoRmax!(Rmax::T, atom::Atom, orbit::Orbit) where T<:Real
    # ==============================================================================
    #  Discretization range in atomic units
    # ==============================================================================
         n = orbit.n
         ℓ = orbit.ℓ
        Zc = atom.Zc
    
        #Rmax = 4(n^2+20)/Zc
        Rmax = Rmax > 0 ? Rmax : 2(2n^2 + 20n + 62)/Zc
        #Rmax = (3.0* n^2 -ℓ*(ℓ+1))/Zc
    
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
    
        Ntot = ID < 3 ? Ntot-1 : I > 3 ? Ntot - 1 : Ntot
        
        h = T(10)/Ntot
        r0 = Rmax / CamiDiff.gridfunction(ID, Ntot, h; p, coords)
    
        return h, r0
    
end

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
    CamiDiff.Grid created: exponential, Float64, Rmax = 16935.0 a.u., Ntot = 3800, h = 0.00263158, r0 = 0.768883

plot_gridfunction(grid, 1:grid.N; title="")
```
The plot is made using CairomMakie.
NB.: `plot_gridfunction` is not part of the `CamiXon` package.
![Image](./assets/exponential_grid.png)
"""
function autoGrid(atom::Atom, orbit::Orbit, T::Type; p=0, coords=[], Ntot=0, Rmax=0, epn=5, k=5, msg=false)

    Rmax = T(Rmax)

    T ∈ [Float64,BigFloat] || println("autoGrid: grid.T = $T => Float64 (was enforced by automatic type promotion)")

    ID = (p < 1) & (length(coords) < 2) ? 1 :
         (p ≥ 1) & (length(coords) < 2) ? (p == 1 ? 3 : 2) :
         (p < 1) & (length(coords) ≥ 2) ? 4 : error("Error: unknown grid")

    Ntot = Ntot == 0 ? autoNtot(orbit) : Ntot
    Rmax = autoRmax!(Rmax, atom, orbit)

    T = T == BigFloat ? T : autoPrecision(Rmax, orbit)
    h, r0 = autoSteps(ID, Ntot, T(Rmax); p, coords)

    return castCamiDiff.Grid(ID, Ntot, T; h, r0, p, coords, epn, k, msg)

end