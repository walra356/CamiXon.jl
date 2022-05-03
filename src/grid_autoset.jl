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
    ID == 3 && return "polynomial"
    ID == 4 && return "linear"

    return error("Error: unknown grid name")

end

# .............. _gridspecs(ID, N, mytype, h, r0; p=5, coords=[0,1]) ...........

function _gridspecs(ID::Int, N::Int, T::Type; h=1, r0=0.001,  p=5, coords=[0,1], epn=5, k=5, msg=true)

    Rmax = ID == 1 ? r0 * _walterjohnson(N, h) :
           ID == 2 ? r0 * _jw_gridfunction(N, h; p) :
           ID == 3 ? r0 * polynomial(coords, h*N)  :
           ID == 4 ? r0 * _linear_gridfunction(N, h) : error("Error: unknown grid type")

    ID = ID ≠ 2 ? ID : p == 1 ? 4 : 2
    name = gridname(ID::Int)
    str_h = repr(h, context=:compact => true)
    str_r0 = repr(r0, context=:compact => true)
    str_Rmax = repr(Rmax, context=:compact => true)
    strA = "Grid created: $(name), $(T), Rmax = "  * str_Rmax * " a.u., Ntot = $N, "

    ID == 1 && return strA * "h = " * str_h * ", r0 = " * str_r0
    ID == 2 && return strA * "p = $p, h = " * str_h * ", r0 = " * str_r0
    ID == 3 && return strA * "coords = $(coords), h = " * str_h * ", r0 = " * str_r0
    ID == 4 && return strA * "p = 1, h = " * str_h * ", r0 = " * str_r0

    return error("Error: unknown grid type")

end

# ....................... autoRmax(atom, orbit) ................................

@doc raw"""
    autoRmax(atom::Atom, orbit::Orbit)

Largest relevant radial distance in a.u. (rule of thumb value)
#### Example:
```
codata = castCodata(2018)
atom = castAtom(Z=1, Q=0, M=1.00782503223, I=1//2, gI=5.585694713; msg=true)
orbit = castOrbit(n=1, ℓ=0)
rmax = autoRmax(atom::Atom, orbit::Orbit); println("rmax = $(rmax) a.u.")
  Atom created: Hydrogen - ¹H (Z = 1, Zc = 1, Q = 0, M = 1.00782503223, I = 1//2, gI = 5.585694713)
  Orbit created: 1s - (n = 1, n′ = 0, ℓ = 0)
  rmax = 63.0 a.u.
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

@doc raw"""
    autoNtot(orbit::Orbit)

Total number of gridpoints (rule of thumb value)
### Example:
```
orbit = castOrbit(1,0)
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

@doc raw"""
    autoPrecision(Rmax::T, orbit::Orbit) where T<:Real

Floating point precision (rule of thumb value)
### Example:
```
atom = castAtom(1)
orbit = castOrbit(1,0)
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

    mytype = Rmax^(ℓ+1) == Inf ? BigFloat : T

    mytype ≠ T && println("\nWarning: overflow protection - changed to BigFloat (Rmax^(ℓ+1) => Inf)\n")

    return mytype

end

# ...................... autoSteps(Ntot, Rmax) .................................

@doc raw"""
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

# ..............................................................................

@doc raw"""
    gridfunction(ID::Int, n::Int, h::T; p=5, coords=[0,1], deriv=0) where T <: Real

* `ID = 1`: exponential grid function,
```math
    f[n] = \text{exp}(h(n-1)) - 1.0
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

# =autoGrid(atom,orbit,codata,T; p=0, coords=[], Nmul=1, epn=7, k=7, msg=true) =

@doc raw"""
    autoGrid(atom::Atom, orbit::Orbit, codata::Codata, T::Type ; p=0, coords=[], Nmul=1, epn=7, k=7, msg=true)

Automatic setting of grid parameters. Important cases: `p=0` (exponential grid
- default), `p=1` (linear grid), `p>1` (quasi-exponential grid)
"""
function autoGrid(atom::Atom, orbit::Orbit, codata::Codata, T::Type ; p=0, coords=[], Nmul=1, epn=7, k=7, msg=true)

    T ∈ [Float64,BigFloat] || println("Warning (autoGrid): grid.T = $T => Float64 (by automatic type promotion)")

    ID = (p < 1) & (length(coords) < 2) ? 1 :
         (p ≥ 1) & (length(coords) < 2) ? (p == 1 ? 4 : 2) :
         (p < 1) & (length(coords) ≥ 2) ? 3 : error("Error: unknown grid")

    Ntot = autoNtot(orbit) * Nmul
    Rmax = autoRmax(atom, orbit)

    T = T == BigFloat ? T : autoPrecision(Rmax, orbit)
    h, r0 = autoSteps(ID, Ntot, Rmax; p, coords)

    return castGrid(ID, Ntot, T; h, r0, p, coords, epn, k, msg)

end
