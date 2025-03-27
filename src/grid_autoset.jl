# SPDX-License-Identifier: MIT

# Copyright (c) 2025 Jook Walraven <69215586+walra356@users.noreply.github.com> and contributors

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# ==============================================================================
#                               grid_autoset.jl
# ==============================================================================

# ....................... autoRmax(atom, orbit) ................................

@doc raw"""
    autoRmax(atom::Atom, n::Int; rmax=0.0)

Largest relevant radial distance in a.u. (rule of thumb value)
```math
    R_{max} = (2.5n^2 + 75.0)/Zc,
```
where ``n`` is the principal quantum number and ``Z_c`` the Rydberg charge
#### Example:
```
julia> atom = castAtom(Z=1, A=1, Q=0);

julia> rmax = autoRmax(atom, n; rmax=0.0); println("rmax = $(rmax) a.u.")
rmax = 77.5 a.u.
```
"""
function autoRmax(atom::Atom, n::Int; rmax=0.0)
# ==============================================================================
#   Discretization range in atomic units
# ==============================================================================
    Zc = atom.Zc
        
    rmax = rmax > 0 ? rmax : (2.5n^2 + 75.0)/Zc
    
    return rmax
        
end

# .......................... autoNtot(n) ...................................

@doc raw"""
    autoNtot(n::Int)

Total number of gridpoints (rule of thumb value)
```math
    N_{tot} = 400 + 100\ n,
```
where ``n`` is the principal quantum number
### Example:
```
julia> autoNtot(1)
500
```
"""
function autoNtot(n::Int)
# ==============================================================================
#   Total number of grid points
# ==============================================================================
    N = 400 + 100 * n
        
    return N
    
end
    
# .................... autoPrecision(rmax, orbit) ..............................

@doc raw"""
    autoPrecision(rmax::T, ℓ = 0) where T<:Real

Floating point precision (rule of thumb value)
### Example:
```
julia> atom = castAtom(Z=1);

julia> orbit = castOrbit(n=1,ℓ=0);

julia> rmax = autoRmax(atom, 0; rmax=0.0); println("rmax = $(rmax) a.u.")
rmax = 77.5 a.u.

julia> o = autoPrecision(rmax, 0); println("precision = $o")
precision = Float64
```
"""
function autoPrecision(rmax::T, ℓ = 0) where T<:Real
# ==============================================================================
#  Floating point precision (rule of thumb value)
# ==============================================================================

    mytype = rmax^(ℓ+1) == Inf ? BigFloat : Float64

    mytype ≠ T && println("\nautoPrecision: rmax^(ℓ+1) => Inf - overflow protection: type promoted to BigFloat\n")

    return mytype

end

@doc raw"""
    autoGrid(atom::Atom, T::Type; p=0, rmax=0, nmax=0, N=0, polynom=[], epn=5, k=5, msg=false) 

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
julia> atom = castAtom(;Z=1, A=1, Q=0, msg=false);

julia> spinorbit = castSpinorbit(n=75, ℓ=0, msg=false);

julia> grid = autoGrid(atom, Float64; nmax=spinorbit.n, msg=true);
Grid: exponential, Float64, rmax = 14137.5, N = 7900, h = 0.00126582, r0 = 0.642684

plot_gridfunction(grid, 1:grid.N; title="")
```
The plot is made using CairomMakie.
NB.: `plot_gridfunction` is not part of the `CamiXon` package.
![Image](../../assets/exponential_grid.png)
"""
function autoGrid(atom::Atom, T::Type; p=0, rmax=0, nmax=0, N=0, polynom=[], epn=5, k=5, msg=false) 

    T ∈ [Float64,BigFloat] || println("autoGrid: grid.T = $T => Float64 (was enforced by automatic type promotion)")
    
    config = collectConfig(atom.config)
    nmax = iszero(nmax) ? maximum([parse(Int, config[i][1]) for i ∈ eachindex(config)]) : nmax
    
    ID = (p < 1) & (length(polynom) < 2) ? 1 :
         (p ≥ 1) & (length(polynom) < 2) ? (p == 1 ? 3 : 2) :
         (p < 1) & (length(polynom) ≥ 2) ? 4 : error("Error: unknown grid")
    
    N = N == 0 ? autoNtot(nmax) : N
    rmax = autoRmax(atom, nmax; rmax)
    
    #T = T == BigFloat ? T : autoPrecision(rmax, orbit)
    
    h = ID == 3 ? T(rmax//N) : T(10//N)
    
    return CamiDiff.castGrid(ID, N, T; h, rmax, p, polynom, epn, k, msg)
    
end
#function autoGrid(atom::Atom, orbit::Orbit, T::Type; p=0, rmax=0, N=0, polynom=[], epn=5, k=5, msg=false)

#    T ∈ [Float64,BigFloat] || println("autoGrid: grid.T = $T => Float64 (was enforced by automatic type promotion)")

#    ID = (p < 1) & (length(polynom) < 2) ? 1 :
#         (p ≥ 1) & (length(polynom) < 2) ? (p == 1 ? 3 : 2) :
#         (p < 1) & (length(polynom) ≥ 2) ? 4 : error("Error: unknown grid") 


#    N = N == 0 ? autoNtot(orbit.n) : N
#    rmax = autoRmax(atom, orbit.n; rmax)

    #T = T == BigFloat ? T : autoPrecision(rmax, orbit)
    
#    h = ID == 3 ? T(rmax//N) : T(10//N)

#    return CamiDiff.castGrid(ID, N, T; h, rmax, p, polynom, epn, k, msg)

# end