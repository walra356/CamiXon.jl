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
    # jw_gridfunction(n, h [; p=5[, deriv=0]]) based on truncated exponential 
    # ==============================================================================
        deriv ≥ 0 || return T(0)
        deriv ≤ p || return T(0)
    
        nul = T(0)
    
        f = deriv > 0 ? h^(deriv)*CamiMath.texp(n*h, nul, p-deriv) : 
                                  CamiMath.texp(n*h, nul, p) - 1  
                                   # note: texp() not exp()
    
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
    
    # ........................ gridname(ID) ........................................
    @doc raw"""
        gridname(ID::Int)
    
    Name corresponding to the grid ID.
    #### Example:
    ```
    n = gridname(2); println("The grid type with ID = 2 is called '$n'.")
      The grid type with ID = 2 is called 'quasi-exponential'.
    ```
    """
    function gridname(ID::Int)
    # ==============================================================================
    #  Name used for `Grid` of given `grid.ID`
    # ==============================================================================
    
        ID == 1 && return "exponential"
        ID == 2 && return "quasi-exponential"
        ID == 3 && return "linear (uniform)"
        ID == 4 && return "polynomial"
    
        return error("Error: unknown grid name")
    
    end
    
    # .............. _gridspecs(ID, N, mytype, h, r0; p=5, coords=[0,1]) ...........
    
    function _gridspecs(ID::Int, N::Int, T::Type; h=1, r0=0.001,  p=5, coords=[0,1], epn=5, k=5, msg=true)
    
        Rmax = ID == 1 ? r0 * _walterjohnson(N, h) :
               ID == 2 ? r0 * _jw_gridfunction(N, h; p) :
               ID == 3 ? r0 * _linear_gridfunction(N, h)  :
               ID == 4 ? r0 * CamiMath.polynomial(coords, h*N) : error("Error: unknown grid type")
    
        ID = ID ≠ 2 ? ID : p == 1 ? 3 : 2
        name = gridname(ID::Int)
        str_h = repr(h, context=:compact => true)
        str_r0 = repr(r0, context=:compact => true)
        str_Rmax = repr(Rmax, context=:compact => true)
        strA = "Grid created: $(name), $(T), Rmax = "  * str_Rmax * " a.u., Ntot = $N, "
    
        ID == 1 && return strA * "h = " * str_h * ", r0 = " * str_r0
        ID == 2 && return strA * "p = $p, h = " * str_h * ", r0 = " * str_r0
        ID == 3 && return strA * "p = 1, h = " * str_h * ", r0 = " * str_r0
        ID == 4 && return strA * "coords = $(coords), h = " * str_h * ", r0 = " * str_r0
    
        return error("Error: unknown grid type")
    
    end
    
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
  [eps(Float64), 0.10517091807564771, 0.22140275816016985, 0.3498588075760032]

grid = castGrid(2, 4, Float64; p = 4, h, r0, msg=true))
grid.r
  create quasi-exponential Grid: Float64, Rmax = 0.491733 a.u., Ntot = 4, p = 4, h = 0.1, r0 = 1.0
  [eps(Float64), 0.10517083333333321, 0.22140000000000004, 0.3498375]

grid = castGrid(3, 4, Float64; coords=[0, 1, 1/2, 1/6, 1/24], h, r0, msg=true)
grid.r
  create polynomial Grid: Float64, Rmax = 0.491733 a.u., Ntot = 4, coords = [0.0, 1.0, 0.5, 0.166666, 0.0416666], h = 0.1, r0 = 1.0
  [eps(Float64), 0.10517083333333334, 0.2214, 0.3498375000000001]

grid = castGrid(4, 4, Float64; h, r0, msg=true)
grid.r
  create linear Grid: Float64, Rmax = 0.4 a.u., Ntot = 4, p = 1, h = 0.1, r0 = 1.0
  [eps(Float64), 0.1, 0.2, 0.3]

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

    r[1] = T == BigFloat ? T("eps(Float64)") : T(eps(Float64))

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
julia> ID = 3; # linear grid
julia> f = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0];
julia> grid = castGrid(ID, length(f), Float64; r0=1.0, h=1.0, k=3, msg=true);
Grid created: linear, Float64, Rmax = 6.0 a.u., Ntot = 6, p = 1, h = 1.0, r0 = 1.0

julia> f′= grid_differentiation(f, grid; k=3); println("f′= $(f′)")
f′= [0.0, 1.9999999999999991, 4.0, 6.000000000000001, 8.0, 10.0]
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
