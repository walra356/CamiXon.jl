# ========== Grid (ID, name, Type, N, r, r′, h, r0, epn, epw, k) ===============

"""
    Grid(ID, name, T, N, r, r′, h, r0, epn, epw, k)

Type with fields:
* `.ID`:   grid identifer name (`::Int`)
* `.name`: grid identifer name (`::String`)
* `.T`:    gridType (`::Type`)
* `.N`:    number of grid points (`::Int`)
* `.r `:   tabulated grid function (`::Vector{T}`)
* `.r′`:   tabulated derivative of grid function (`::Vector{T}`)
* `.h` :   grid step multiplyer (`::T`)
* `.r0`:   grid scale factor (`::T`)
* `.epn`:  number of endpoints used for trapezoidal endpoint correction (must be odd) (`::Int`)
* `.epw`:  trapezoidal endpoint weights for n=1:epn (`::Vector{Vector{T}}`)
* `.k`:    Adams-Moulton order (`::Int`)
The object `Grid` is best created with the function [`castGrid`](@ref).
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

# ====== castGrid(ID, T, N; h=1, r0=0.01,  p=5, coords=[0,1], epn=5, k=7) ====

"""
    castGrid(ID::Int, N::Int, T::Type; h=1, r0=1,  p=5, coords=[0,1], epn=5, k=7, msg=false)

Method to create the Grid object

`ID = 1`: exponential grid,
`ID = 2`: quasi-exponential grid,
`ID = 3`: linear grid
`ID = 4`: polynomial grid
#### Examples:
```
h = 0.1
r0 = 1.0
grid = castGrid(1, 4, Float64; h, r0, msg=true)
grid.r
  create exponential Grid: Float64, Rmax = 0.491825 a.u., Ntot = 4, h = 0.1, r0 = 1.0
  [1.0e-100, 0.10517091807564771, 0.22140275816016985, 0.3498588075760032]

grid = castGrid(2, 4, Float64; p = 4, h, r0, msg=true))
grid.r
  create quasi-exponential Grid: Float64, Rmax = 0.491733 a.u., Ntot = 4, p = 4, h = 0.1, r0 = 1.0
  [1.0e-100, 0.10517083333333321, 0.22140000000000004, 0.3498375]

grid = castGrid(3, 4, Float64; coords=[0, 1, 1/2, 1/6, 1/24], h, r0, msg=true)
grid.r
  create polynomial Grid: Float64, Rmax = 0.491733 a.u., Ntot = 4, coords = [0.0, 1.0, 0.5, 0.166666, 0.0416666], h = 0.1, r0 = 1.0
  [1.0e-100, 0.10517083333333334, 0.2214, 0.3498375000000001]

grid = castGrid(4, 4, Float64; h, r0, msg=true)
grid.r
  create linear Grid: Float64, Rmax = 0.4 a.u., Ntot = 4, p = 1, h = 0.1, r0 = 1.0
  [1.0e-100, 0.1, 0.2, 0.3]

grid.r′
  [0.1, 0.1, 0.1, 0.1]
```
"""
function castGrid(ID::Int, N::Int, T::Type; h=1, r0=0.001,  p=5, coords=[0,1], epn=5, k=5, msg=false)
# ==============================================================================
#  castGrid: creates the grid object
# ==============================================================================
    h = convert(T, h)
    r0 = convert(T, r0)
    coords = convert.(T, coords)
    epw = [convert.(T, trapezoidal_epw(n; rationalize=true)) for n=1:2:epn]
    name = gridname(ID)

    r = r0 * [gridfunction(ID, n-1, h; p, coords) for n=1:N]
    r′= r0 * [gridfunction(ID, n-1, h; p, coords, deriv=1) for n=1:N]     # r′= dr/dn

    r[1] = T == BigFloat ? T("1.0e-100") : T(1.0e-100)

    msg && println(_gridspecs(ID, N, T; h, r0,  p, coords, epn, k, msg))

    return Grid(ID, name, T, N, r, r′, h, r0, epn, epw, k)

end

# =============== findIndex(rval, grid) ========================================

@doc raw"""
    findIndex(rval::T, grid::Grid{T}) where T<:Number

The grid index corresponding to the position `rval` on the `grid`.
#### Example:
```
h = 0.1
r0 = 1.0
grid = castGrid(1, 4, Float64; h, r0)
r = grid.r; println("r[3] = $(r[3])")
  Grid created: exponential, Float64, Rmax = 0.491825 a.u., Ntot = 4, h = 0.1, r0 = 1.0
  r[3] = 0.22140275816016985

findIndex(0.222, grid)
  3
```
"""
function findIndex(rval::T, grid::Grid{T}) where T<:Number
# ==============================================================================
#  grid index of rval, e.g., rval -> classical turning point
# ==============================================================================
    N = grid.N
    r = grid.r

    r[1] ≤ rval ≤ r[end] || error("Error: radial distance in a.u. outside grid range")

    n = N
    while rval < r[n]     # below classical threshhold
        n > 1 ? n -= 1 : break
    end

    return n

end

# =============== grid_differentiation(f, grid; k=3)) ======================

@doc raw"""
    grid_differentiation(f::Vector{T}, grid::Grid{T}; k=3) where T<:Real

``k^{th}``-order lagrangian *differentiation* of the analytic function ``f``,
tabulated in forward order on a [`Grid`](@ref) of ``n`` points, ``f[1:n]``.
#### Example:
```
ID = 4 # linear grid
f = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0]
grid = castGrid(ID, length(f), Float64; r0=1.0, h=1.0, k=3)  # linear grid
f′= grid_differentiation(f, grid; k=3); println("f′= $(f′)")
  Grid created: linear, Float64, Rmax = 6.0 a.u., Ntot = 6, p = 1, h = 1.0, r0 = 1.0
  f′= [0.0, 2.0, 4.0, 6.0, 7.999999999999998, 9.999999999999993]
```
"""
function grid_differentiation(f::Vector{T}, grid::Grid{T}; k=3) where T<:Real

    r′= grid.r′

    l = length(f)
    f′= [fdiff_differentiation(f, T(v); k) for v=1:l]

    return f′ ./ r′

end

# =============== grid_integration(f, n1, n2, grid) ===================

@doc raw"""
    grid_integration(f::Vector{T}, n1::Int, n2::Int, grid::Grid{V}) where {T<:Real, V<:Real}

Integral of the function ``f=[f_0,⋯\ f_n]`` tabulated on a [`Grid`](@ref)
using the trapezoidal rule optimized with endpoint correction by the
weightsvector `grid.epw`,
```math
    ∫_{0}^{r_n} f(r) dr = ∫_{0}^{n} f(x) r^{\prime}(x) dx,
```
where the latter integral corresponds to the optimized trapezoidal rule for a
uniform grid (see [`trapezoidal_integration`](@ref)). The rule is exact for
polynonials of degree ``d=0,\ 1,⋯\ k-1``, where ``k=`` `grid.epn`.
For ``k=1`` the rule reduces to the ordinary trapezoidal rule (weights = [1/2]).
#### Example:
```
f1s(r) = 2.0*r*exp(-r);  # hydrogen 1s wavefunction (reduced and unit normalized)
N = 1000;
grid = castGrid(1, N, Float64; h=0.01, r0=0.005)
    create exponential Grid: Float64, Rmax = 110.127 (a.u.), Ntot = 1000, h = 0.01, r0 = 0.005

r = grid.r;
f2 = [f1s(r[n])^2 for n=1:N];
grid_integration(f2, 1:N, grid) == grid_integration(f2, 1, N, grid)
    true

norm = grid_integration(f2, 1:N, grid)

    1.0
```
"""
function grid_integration(f::Vector{T}, n1::Int, n2::Int, grid::Grid{V}) where {T<:Real, V<:Real}
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
        epn = Base.isodd(n÷2) ? (n÷2) : n÷2-1           # endpoint number
        epi = Base.isodd(n÷2) ? (n÷2)÷2+1 : (n÷2-1)÷2+1 # index endpoint weights
    end

    w = Base.ones(T,n)
    w[1:epn] = epw[epi]
    w[end-epn+1:end] = Base.reverse(epw[epi])

    return LinearAlgebra.dot(f .* r′, w)

end
function grid_integration(f::Vector{T}, itr::UnitRange, grid::Grid{V}) where {T<:Real, V<:Real}

    return grid_integration(f, itr.start, itr.stop, grid)

end
