# ======================== Grid (ID, name, Type, N, r, r′, h, r0, epn, epw, k) ===============

"""
    Grid(ID, name, T, N, r, r′, h, r0, epn, epw, k)

Type with fields:
* `      .ID`::Int                 grid identifer name
* `    .name`::String              grid identifer name
* `       .T`::Type                gridType
* `       .N`::Int                 number of grid points
* `      .r `::Vector{T}           tabulated grid function
* `      .r′`::Vector{T}           tabulated derivative of grid function
* `      .h` ::T                   grid step multiplyer
* `      .r0`::T                   grid scale factor
* `     .epn`::Int                 number of endpoints used for trapezoidal endpoint correction (must be odd)
* `     .epw`::Vector{Vector{T}}   trapezoidal endpoint weights for n=1:epn
* `       .k`::Int                 Adams-Moulton order

The type `Grid` is best created by the function `createGrid`.
"""
struct Grid{T}
    ID::Int
    name::String
    T::Type
    N::Int
    r ::Vector{T}
    r′::Vector{T}
    h::T
    r0::T
    epn::Int
    epw::Vector{Vector{T}}
    k::Int
end

# ======================== gridfunction(n, h; deriv=0) =========================

function _walterjohnson(n::Int, h::T; deriv=0) where T <: Real
# ==============================================================================
#  gridfunction(n, h) = (exp((n-1) * h)-1.0) # gridfunction from Walter Johnson
# ==============================================================================
    deriv ≥ 0 || return 0.0

    f = deriv > 0 ? h^(deriv)*exp(n*h) : exp(n*h)-1

    return f

end

# ..............................................................................

function _jw_gridfunction(n::Int, h::T; p=5, deriv=0) where T <: Real
# ==============================================================================
# jw_gridfunction(n, h [; p=5[, deriv=0]]) based on truncated exponential texp()
# ==============================================================================
    deriv ≥ 0 || return T(0)
    deriv ≤ p || return T(0)

    nul = T(0)

    f = deriv > 0 ? h^(deriv)*texp(n*h, nul, p-deriv) : texp(n*h, nul, p) - 1  # note: texp() not exp()

    return f

end

# ..............................................................................

function _linear_gridfunction(n::Int, h::T; deriv=0) where T <: Real
# ==============================================================================
#  linear_gridfunction(n, h; deriv) = n * h
# ==============================================================================
    deriv ≥ 0 || return T(0)
    deriv ≤ 1 || return T(0)

    f = deriv > 0 ? h : h * n

    return f

end

# ..............................................................................

@doc raw"""
    gridfunction(ID::Int, n::Int, h::T [; p=5, coords=[0,1], deriv=0]) where T <: Real

* `ID = 1`: exponential grid function,
```math
    f[n] = exp(h(n-1)) - 1.0
```
* `ID = 2`: quasi-exponential grid function (linear grid for p = 1),
```math
    f[n] = h(n-1) + \frac{1}{2}(h(n-1))^2 + \cdots + \frac{1}{p!}(h(n-1))^p
```
* `ID = 3`: polynomial grid function based on `polynom` ``c = [c_1,c_2,\ldots,c_p]``,
```math
    f[n] = c_1h(n-1) + c_2(h(n-1))^2 + \cdots + c_p(h(n-1))^p
```
* `ID = 4`: linear grid function,
```math
    f[n] = (n-1) * h
```
#### Examples:
```
h = 0.1
r = [gridfunction(1, n-1, h) for n=1:5]                            # exponential
 [0.0, 0.10517091807564771, 0.22140275816016985, 0.3498588075760032, 0.49182469764127035]

r = [gridfunction(2, n-1, h; p = 4) for n=1:5]  # quasi exponential (degree p=4)
 [0.0, 0.10517083333333321, 0.22140000000000004, 0.3498375, 0.49173333333333336]

r = [gridfunction(3, n-1, h; coords = [0,1,1/2,1/6,1/24]) for n=1:5]  # polynomial (degree p=4)
 [0.0, 0.10517083333333334, 0.2214, 0.3498375000000001, 0.49173333333333336]

r = [gridfunction(4, n-1, h) for n=1:5]              # linear
  [0.0, 0.1, 0.2, 0.3, 0.4]

r′= [gridfunction(4, n-1, h; deriv=1) for n=1:5]     # linear (first derivative)
   [0.1, 0.1, 0.1, 0.1, 0.1]
```
"""
function gridfunction(ID::Int, n::Int, h::T; p=5, coords=[0,1], deriv=0) where T <: Real

    ID == 1 && return _walterjohnson(n, h; deriv)
    ID == 2 && return _jw_gridfunction(n, h; deriv, p)
    ID == 3 && return polynomial(coords, h*n; deriv)
    ID == 4 && return _linear_gridfunction(n, h; deriv)

    return error("Error: unknown gridfunction")

end

# ====== createGrid(ID, T, N; h=1, r0=0.01,  p=5, coords=[0,1], epn=7, k=7) ====

@doc raw"""
    createGrid(ID::Int, N::Int, T::Type [; h=1, r0=0.01,  p=5, coords=[0,1], epn=7, k=7, msg=true]))

Create the Grid object

`ID = 1`: exponential grid,
`ID = 2`: quasi-exponential grid,
`ID = 3`: polynomial grid
`ID = 4`: linear grid
#### Examples:
```
h = 0.1
r0 = 1.0
grid = createGrid(1, 4, Float64; h, r0)                 # exponential grid
grid.r
 [0.0, 0.10517091807564771, 0.22140275816016985, 0.3498588075760032]

grid = createGrid(2, 4, Float64; p = 4, h, r0)          # quasi-exponential grid
grid.r
 [0.0, 0.10517083333333321, 0.22140000000000004, 0.3498375]

grid = createGrid(3, 4, Float64; coords=[0, 1, 1/2, 1/6, 1/24], h, r0)  # polynomial grid
grid.r
 [0.0, 0.10517083333333334, 0.2214, 0.3498375000000001]

grid = createGrid(4, 4, Float64; h, r0)                 # linear grid
grid.r
 [0.0, 0.1, 0.2, 0.3]
grid.r′
 [0.1, 0.1, 0.1, 0.1]
```
"""
function createGrid(ID::Int, N::Int, T::Type; h=1, r0=0.001,  p=5, coords=[0,1], epn=7, k=7, msg=true)
# ==============================================================================
#  createGrid: creates the grid object
# ==============================================================================
    h = myconvert(T, h)
    r0 = myconvert(T, r0)
    coords = myconvert.(T, coords)
    epw = [myconvert.(T, trapezoidal_weights(n; rationalize=true)) for n=1:2:epn]
    name = gridname(ID)

    r = r0 * [gridfunction(ID, n-1, h; p=p, coords) for n=1:N]
    r′= r0 * [gridfunction(ID, n-1, h; p=p, coords, deriv=1) for n=1:N]

    msg && println(gridspecs(ID, N, T, h, r0; p, coords))

    return Grid(ID, name, T, N, r, r′, h, r0, epn, epw, k)

end

# ....................... autoRmax(atom, orbit) ................................

"""
    autoRmax(atom::Atom, orbit::Orbit)

Discretization range in atomic units (rule of thumb value)
### Example:
```
atom = createAtom(1)
orbit = createOrbit(1,0)
autoRmax(atom, orbit)
 Atom created: Hydrogen - ¹H (Z = 1, Zc = 1, Q = 0, M = 1.0, I = 1//2, gI = 5.5)
 Orbit created: 1s - (n = 1, n′ = 0, ℓ = 0)

 63.0
```
"""
function autoRmax(atom::Atom, orbit::Orbit)
# ==============================================================================
#  Discretization range in atomic units
# ==============================================================================
     n = orbit.n
    Zc = atom.Zc

    Rmax = 3(n^2+20)/Zc

    return Rmax

end

# .......................... autoNtot(orbit) ...................................

"""
    autoNtot(orbit::Orbit)

Total number of gridpoints (rule of thumb value)
### Example:
```
orbit = createOrbit(1,0)
autoNtot(orbit)
 Orbit created: 1s - (n = 1, n′ = 0, ℓ = 0)

 100
```
"""
function autoNtot(orbit::Orbit)
# ==============================================================================
#  Total number of grid points
# ==============================================================================

    n = orbit.n

    Ntot = (50 + 50 * n)

    return Ntot

end

# .................... autoPrecision(Rmax, orbit) ..............................

"""
    autoPrecision(Rmax::T, orbit::Orbit) where T<:Real

Floating point precision (rule of thumb value)
### Example:
```
atom = createAtom(1)
orbit = createOrbit(1,0)
Rmax = autoRmax(atom, orbit)
autoPrecision(Rmax, orbit)
 Atom created: Hydrogen - ¹H (Z = 1, Zc = 1, Q = 0, M = 1.0, I = 1//2, gI = 5.5)
 Orbit created: 1s - (n = 1, n′ = 0, ℓ = 0)

 Float64
```
"""
function autoPrecision(Rmax::T, orbit::Orbit) where T<:Real
# ==============================================================================
#  Floating point precision (rule of thumb value)
# ==============================================================================

    ℓ = orbit.ℓ

    mytype = Rmax^(ℓ+1) == Inf ? BigFloat : Float64

    return mytype

end

# ...................... autoSteps(Ntot, Rmax) .................................

"""
    autoSteps(ID::Int, Ntot::Int, Rmax::T; p=5, coords=[0,1]) where T<:Real

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


# =============== grid_trapezoidal_integral(f, n1, n2, grid) ===================

"""
    grid_trapezoidal_integral(f::Vector{T}, n1::Int, n2::Int, grid::Grid{T}) where T<:Real

Generalized trapezoidal integral with endpoint correction on `epn = grid.epn points.
#### Example:
```
f1s(r) = 2.0*r*exp(-r)  # hydrogen 1s wavefunction (reduced and unit normalized)
N = 1000
grid = createGrid(1, N, Float64; h=0.01, r0=0.005)
r = grid.r
f2 = [f1s(r[n])^2 for n=1:N]
norm = grid_trapezoidal_integral(f2, 1:N, grid)
 create exponential Grid: Float64, Rmax = 110.127 (a.u.), Ntot = 1000, h = 0.01, r0 = 0.005

 1.0
```
"""
function grid_trapezoidal_integral(f::Vector{T}, n1::Int, n2::Int, grid::Grid{T}) where T<:Real
# ==============================================================================
#  trapezoidal integral over the grid indices [n1:n2] with 1 ≤ n1,n2 ≤ N
# ==============================================================================
    f = f[n1:n2]
    r′= grid.r′[n1:n2]
    n = n2 - n1 + 1

    n > 1 || return 0

    epn = grid.epn   # endpoint number
    epw = grid.epw   # endpoint weights array

    if n ≥ 2epn
        epi = epn÷2+1        # index endpoint weights
    else
        epn = Base.isodd(n÷2) ? (n÷2) : n÷2-1              # endpoint number
        epi = Base.isodd(n÷2) ? (n÷2)÷2+1 : (n÷2-1)÷2+1    # index endpoint weights
    end

    w = Base.ones(T,n)
    w[1:epn] = epw[epi]
    w[end-epn+1:end] = Base.reverse(epw[epi])

    return LinearAlgebra.dot(f .* r′, w)

end
function grid_trapezoidal_integral(f::Vector{T}, itr::UnitRange, grid::Grid{T}) where T<:Real

    return grid_trapezoidal_integral(f, itr.start, itr.stop, grid)

end
