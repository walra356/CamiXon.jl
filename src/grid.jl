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
The object `Grid` is best created by the function [`castGrid`](@ref).
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

# ====== castGrid(ID, T, N; h=1, r0=0.01,  p=5, coords=[0,1], epn=7, k=7) ====

"""
    castGrid(ID::Int, N::Int, T::Type; h=1, r0=1,  p=5, coords=[0,1], epn=7, k=7, msg=true)

Method to create the Grid object

`ID = 1`: exponential grid,
`ID = 2`: quasi-exponential grid,
`ID = 3`: polynomial grid
`ID = 4`: linear grid
#### Examples:
```
h = 0.1
r0 = 1.0
grid = castGrid(1, 4, Float64; h, r0)
grid.r
  create exponential Grid: Float64, Rmax = 0.491825 a.u., Ntot = 4, h = 0.1, r0 = 1.0
  [0.0, 0.10517091807564771, 0.22140275816016985, 0.3498588075760032]

grid = castGrid(2, 4, Float64; p = 4, h, r0)
grid.r
  create quasi-exponential Grid: Float64, Rmax = 0.491733 a.u., Ntot = 4, p = 4, h = 0.1, r0 = 1.0
  [0.0, 0.10517083333333321, 0.22140000000000004, 0.3498375]

grid = castGrid(3, 4, Float64; coords=[0, 1, 1/2, 1/6, 1/24], h, r0)
grid.r
  create polynomial Grid: Float64, Rmax = 0.491733 a.u., Ntot = 4, coords = [0.0, 1.0, 0.5, 0.166666, 0.0416666], h = 0.1, r0 = 1.0
  [0.0, 0.10517083333333334, 0.2214, 0.3498375000000001]

grid = castGrid(4, 4, Float64; h, r0)
grid.r
  create linear Grid: Float64, Rmax = 0.4 a.u., Ntot = 4, p = 1, h = 0.1, r0 = 1.0
  [0.0, 0.1, 0.2, 0.3]

grid.r′
  [0.1, 0.1, 0.1, 0.1]
```
"""
function castGrid(ID::Int, N::Int, T::Type; h=1, r0=0.001,  p=5, coords=[0,1], epn=5, k=5, msg=true)
# ==============================================================================
#  castGrid: creates the grid object
# ==============================================================================
    h = convert(T, h)
    r0 = convert(T, r0)
    coords = convert.(T, coords)
    epw = [convert.(T, trapezoidal_weights(n; rationalize=true)) for n=1:2:epn]
    name = gridname(ID)

    r = r0 * [gridfunction(ID, n-1, h; p, coords) for n=1:N]
    r′= r0 * [gridfunction(ID, n-1, h; p, coords, deriv=1) for n=1:N]

    msg && println(_gridspecs(ID, N, T; h, r0,  p, coords, epn, k, msg))

    return Grid(ID, name, T, N, r, r′, h, r0, epn, epw, k)

end

# =============== grid_lagrange_derivative(f, grid; k=5)) ====

@doc raw"""
    grid_lagrange_derivative(f::Vector{T}, grid::Grid{T}; k=5) where T<:Real

``k^{th}``-order lagrangian *differentiation* of the analytic function ``f``,
tabulated in forward order on a [`Grid`](@ref) of ``n`` points, ``f[1],\ ⋯,
\ f[n]``; ``m`` is the multiplier for intermediate positions (for ``m=1``
*without* intermediate points).
#### Example:
```
ID = 4 # linear grid
f = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0]
grid = castGrid(ID, length(f), Float64; r0=1.0, h=1.0, k=3)  # linear grid
f′= grid_lagrange_derivative(f, grid, k=4)
f′= ceil.(f′;sigdigits=2); println(f′)
  create linear Grid: Float64, Rmax = 11.0 (a.u.), Ntot = 11, p = 1, h = 1.0, r0 = 1.0
  [0.0, 2.0, 4.0, 6.0, 8.1, 11.0, 12.0, 14.0, 17.0, 18.0, 20.0]
```
"""
function grid_lagrange_derivative(f::Vector{T}, grid::Grid{T}; k=5) where T<:Real

    N = grid.N
    r′= grid.r′

    ∇ = CamiXon.f_diff_weights_array(k)
    β = [CamiXon.f_diff_expansion_coeffs_differentiation(k, x) for x=-k:0]
    w = [CamiXon.bwd_diff_expansion_weights(β[i], ∇) for i ∈ Base.eachindex(β)]
    u = Base.append!(repeat(w[1:1],N-k-1),w)
    v = CamiXon.f_diff_function_sequences(f , k, 1)

    f′= [u[i] ⋅ v[i] for i ∈ Base.eachindex(u)]

    return f′ ./ r′

end

# =============== grid_trapezoidal_integral(f, n1, n2, grid) ===================

@doc raw"""
    grid_trapezoidal_integral(f::Vector{T}, n1::Int, n2::Int, grid::Grid{T}) where T<:Real

Integral of the function ``f=[f_0,⋯,\ f_n]`` tabulated on a [`Grid`](@ref)
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
f1s(r) = 2.0*r*exp(-r)  # hydrogen 1s wavefunction (reduced and unit normalized)
N = 1000
grid = castGrid(1, N, Float64; h=0.01, r0=0.005)
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
