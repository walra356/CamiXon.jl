include atom.jl
# ================== Pos(Na, Nlctp, Nmin, Nuctp, Nb, N, nodes) =================
"""
    Pos(Na::Int, Nlctp::Int, Nmin::Int, Nuctp::Int, Nb::Int, N::Int, nodes::Int)

Type with fields:
* `   .Na`: grid index of last leading point (`::Int`)
* `.Nlctp`: grid index of classical turning point (`::Int`)
* ` .Nmin`: grid index of (screened) potential minimum (`::Int`)
* `.Nuctp`: grid index of classical turning point (`::Int`)
* `   .Nb`: grid index first trailing point (`::Int`)
* `    .N`: grid index last point (`::Int`)
* `.nodes`: number of nodes  (`::Int`)

Mutable struct to hold special grid indices as well as the number of nodes
#### Examples:
```
pos = Pos(1, 2, 3, 4, 5, 6, 7)
pos.Nuctp
 4

pos.Nuctp = 8
pos
 Pos(1, 2, 3, 8, 5, 6, 7)
```
"""
mutable struct Pos
    Na::Int
    Nlctp::Int
    Nmin::Int
    Nuctp::Int
    Nb::Int
    N::Int
    nodes::Int
end

# ===========   Grid (ID, name, Type, N, r, r′, h, r0, epn, epw, k) ============

"""
    Def(T, atom, orbit, pot, scr, o1, o2, o3, pos, epn, k, am, matLD)

Type with fields:
* `    .T`: gridType (`::Type`)
* ` .atom`: atom object (`::Atom`)
* `.orbit`: orbit object (`::Orbit`)
* `  .pot`: tabulated potential function (`::Vector{T}`)
* `  .scr`: tabulated screening function (`::Vector{T}`)
* `   .o1`: vector of zero-filled matrices (`::Vector{Matrix{T}}`)
* `   .o2`: vector of zero-filled matrices (`::Vector{Matrix{T}}`)
* `   .o3`: vector of unit-filled matrices (`::Vector{Matrix{T}}`)
* `  .pos`: object containing Na, Nlctp, Nmin, Nuctp, Nb, N and nodes (`::Pos`)
* `  .epn`: number of endpoints trapezoidal correction - must be odd (`::Int`)
* `    .k`: Adams-Moulton order (`::Int`)
* `   .am`: Adams-Moulton weight coefficients (`::Vector{T}`)
* `.matLD`: Lagrangian differentiation matrix (`::Matrix{T} )

The object `Def` is best created by the function `castDef`.
"""
struct Def{T}
    T::Type
    atom::Atom
    orbit::Orbit
    pot::Vector{T}          # tabulated potential function
    scr::Vector{T}          # tabulated screening function
    o1::Vector{Matrix{T}}   # vector of zero-filled matrices
    o2::Vector{Matrix{T}}   # vector of zero-filled matrices
    o3::Vector{Matrix{T}}   # vector of unit-filled matrices
    pos::Pos                # object containing Nmin, Na, Nuctp, Nb, N and nodes
    epn::Int                # number of endpoints trapezoidal correction
    k::Int                  # Adams-Moulton order
    am::Vector{T}           # Adams-Moulton weight coefficients
    matLD::Matrix{T}        # Lagrangian differentiation matrix
end

"""
    castDef(grid::Grid{T}, atom::Atom, orbit::Orbit) where T <: Real

Create the Def object starting from the Grid and atomic properties.
"""
function castDef(grid::Grid{T}, atom::Atom, orbit::Orbit; scr=nothing) where T <: Real
# ================================================================================
# castDef(grid, atom, orbit) # reference arrays
# ================================================================================
    N = grid.N
    r = grid.r
    k = grid.k
    epn = grid.epn
    Z = atom.Z
    ℓ = orbit.ℓ

    r[N]^(ℓ+1) < Inf || error("Error: numerical overflow (Inf)")

    Z = myconvert(T, Z)
    num = myconvert(T, ℓ*(ℓ + 1)//2)

    r1 = T(1.0e-100)  # quasi zero
    pot = ℓ > 0 ? [(-Z + num/r[n])/r[n] for n=1:N] : [-Z/r[n] for n=1:N]
    pot[1] = ℓ > 0 ? (-Z + num/r1)/r1 : -Z/r1
    pot = convert.(T,pot)
    scr = isnothing(scr) ? zeros(T,N) : scr
    o1 = [fill(myconvert(T,0), (2,2)) for n=1:N]
    o2 = [fill(myconvert(T,0), (2,2)) for n=1:N]
    o3 = [fill(myconvert(T,1), (2,2)) for n=1:N]
    pos = Pos(k+1, 0, 1, 0, N-k, N, 0)  # Pos(Na, Nlctp, Nmin, Nuctp, Nb, N, nodes)
    am = myconvert.(T, create_adams_moulton_weights(k; rationalize=true))
    matLD = myconvert.(T, create_lagrange_differentiation_matrix(k))

    return Def(T, atom, orbit, pot, scr, o1, o2, o3, pos, epn, k, am, matLD)

end

@doc raw"""
    initE(def::Def{T}; E=nothing) where T<:Real

Autogenerated seed value for the energy (option: E as a manual seed)
#### Example:
```
codata = castCodata(2018)
atom = castAtom(Z=1, Q=0, M=1.00782503223, I=1//2, gI=5.585694713)
orbit = castOrbit(n=1, ℓ=0)
Ecal = convert(Float64, bohrformula(atom.Z, orbit.n))
grid = autoGrid(atom, orbit, codata; msg=false)
def = castDef(grid, atom, orbit)
  Atom created: Hydrogen - ¹H (Z = 1, Zc = 1, Q = 0, M = 1.00782503223, I = 1//2, gI = 5.585694713)
  Orbit created: 1s - (n = 1, n′ = 0, ℓ = 0)

E = initE(def); println("E = $E")
  E = -0.03508495857961283

E = initE(def; E=Ecal); println("E = $E")
  E = -0.5
```
"""
function initE(def::Def{T}; E=nothing) where T<:Real

    isnothing(E) || return myconvert(T, E)

    N = def.pos.N
    ℓ = def.orbit.ℓ
    v = def.pot
    s = def.scr

    pot = v .+ s

    Emax = pot[N]
    Emin = minimum(pot[2:N])

    E = iszero(ℓ) ? 2.0Emax : 0.9Emin

    return E

end

function _get_Na(Z::Vector{Complex{T}}, def::Def{T}) where T<:Real
# ==============================================================================
#  grid index of starting point outward numerical integration
# ==============================================================================
    k = def.k

    ref = T(1.0e-10)

    Na = findfirst(x -> abs(x) > ref, real(Z))
    Na = isnothing(Na) ? k+1 : Na > 0 ? max(k+1, Na) : k+1

    return Na

end

function _get_Nb(Z::Vector{Complex{T}}, def::Def{T}) where T<:Real
# ==============================================================================
#  grid index of starting point inward numerical integration
# ==============================================================================
    k = def.k
    N = def.pos.N

    ref = T(1.0e-10)

    Nb = findlast(x -> abs(x) > ref, real(Z))
    Nb = isnothing(Nb) ? N-k : Nb > 0 ? min(N-k, Nb) : N-k

    return Nb

end

function _get_Nmin(def::Def{T}) where T<:Real
# ==============================================================================
#  grid index of potential minimum
# ==============================================================================
    N = def.pos.N
    v = def.pot
    s = def.scr

    pot = v .+ s

    n = 1
    while pot[n+1] < pot[n]
        n < N-1 ? n += 1 : break
    end

    return n

end

function _get_Nuctp(E::T, def::Def{T}) where T<:Real
# ==============================================================================
#  grid index of upper classical turning point
# ==============================================================================
    Nmin = def.pos.Nmin
    N = def.pos.N
    v = def.pot
    s = def.scr

    pot = v .+ s

    n = N-1

    E < pot[Nmin] && println("Warning: E < Emin - force Nuctp = Nmin")

    while E < pot[n]     # below classical threshhold
        n > Nmin ? n -= 1 : break
    end

    n < N || error("Error: (outer) classical turning point outside grid")

    return n

end

function _get_Nlctp(E::T, def::Def{T}) where T<:Real
# ==============================================================================
#  grid index of lower classical turning point
# ==============================================================================
    Nmin = def.pos.Nmin
    N = def.pos.N
    v = def.pot
    s = def.scr
    ℓ = def.orbit.ℓ

    ℓ > 0 || return 0

    pot = v .+ s

    E < pot[Nmin] && println("Warning: E < Emin - force Nuctp = Nmin")

    n = 2

    while E < pot[n]     # below classical threshhold
        n < Nmin ? n += 1 : break
    end

    n < N || error("Error: inner classical turning point outside grid")

    return n

end

function matG(E::T, grid::Grid{T}, def::Def{T}) where T<:Real
# ==============================================================================
# matG - coupling matrix - Johnson (2.54)
# ==============================================================================
    r′= grid.r′
    o = def.o1
    v = def.pot
    s = def.scr

    pot = v .+ s

    two = T(2)

    for n ∈ eachindex(o)
        o[n][1,2] =  r′[n]
        o[n][2,1] =  r′[n] * two * (-E + pot[n])
    end

    return o

end

function matσ(E::T, grid::Grid{T}, def::Def{T}) where T<:Real
# ==============================================================================
# matσ - coupling matrix - Johnson (2.54)
# ==============================================================================
    r = grid.r
    r′= grid.r′
    o = def.o2
    Z = def.atom.Z
    ℓ = def.orbit.ℓ
    s = def.scr

    Zet = T(Z)
    one = T(1)
    two = T(2)
    num = -T(ℓ + 1)

    for n ∈ eachindex(o)
        o[n][1,2] = r′[n]                                  # b in Johnson (2.69)
        o[n][2,1] = r′[n] * two * (-E - Zet/r[n] + s[n])   # c in Johnson (2.69)
        o[n][2,2] = r′[n] * two * num / r[n]               # d in Johnson (2.69)
    end

    r1 = T(1.0e-100)  # quasi zero
    o[1][2,1] = r′[1] * two * (-E - Zet / r1 + s[1])
    o[1][2,2] = r′[1] * two * num / r1

    return o

end

function matMinv(E::T, grid::Grid{T}, def::Def{T}, amEnd::T) where T<:Real
# ==============================================================================
# matMinv - Adams-Moulton correction matrix - Johnson (2.56)
# ==============================================================================
    N = grid.N
    r′= grid.r′
    v = def.pot
    s = def.scr
    o = def.o3

    one = T(1)
    two = T(2)

    λ::Vector{T} = r′ * amEnd       # note that am_0=c[end] (the weights are stored in reverse order)

    for n ∈ eachindex(o)
        o[n][1,2] = λ[n]
        o[n][2,1] = λ[n]  * two * (-E + v[n] + s[n])
    end

    Δ = [one - o[n][1,2]*o[n][2,1] for n ∈ eachindex(o)]

    return o ./ Δ

end

# ============================= OUTSCH sector ==================================

function _matOUTSCH(k::Int, matLD::Matrix{T}, matσ::Vector{Matrix{T}}) where T<:Real

    mat = fill(T(0),(2k,2k))

    for i=1:k
        for j=1:k
            mat[i,j] = matLD[1+i,1+j]
            mat[k+i,k+j] = mat[i,j]
        end
    end

    for i=1:k  mat[i,k+i] = -matσ[i][1,2]  end
    for i=1:k  mat[k+i,i] = -matσ[i][2,1]  end
    for i=1:k  mat[k+i,k+i] = mat[k+i,k+i] - matσ[i][2,2]  end

    return mat

end
# ..................................................................................
function _vecOUTSCH(k::Int, matLD::Matrix{T}, p::T, q::T) where T<:Real

        v = [matLD[2,1]]

    for i=2:k  push!(v, matLD[1+i,1] * p)  end       # Johnson (2.67)
    for i=1:k  push!(v, matLD[1+i,1] * q)  end       # Johnson (2.68)

    return -v

end
# ..................................................................................
function _update_Z!(Z::Vector{Complex{T}}, o::Vector{Complex{T}}, Na::Int, k::Int) where T<:Real

    for n=Na:Na+k
        Z[n] = o[n-Na+1]
    end

    return Z

end
# ..................................................................................
function OUTSCH(grid::Grid{T}, def::Def{T}, σ::Vector{Matrix{T}}) where T<:Real

        r = grid.r
        k = def.k
        N = def.pos.N
        ℓ = def.orbit.ℓ
     Zval = def.atom.Z
    matLD = def.matLD

    p = T(1)
    q = myconvert(T, -Zval//(ℓ + 1) )
    num = T(ℓ + 1)

    o = zeros(Complex{T}, k+1)
    o[1] = p + q * im                       # boundary condition

    m = _matOUTSCH(k, matLD, σ)
    v = _vecOUTSCH(k, matLD, p, q)

    u = inv(m) * v                          # solve 2dx2d set of equations

    for n=1:k
        o[n+1] = u[n] + u[n+k]*im
    end

    Z = zeros(Complex{T},N)

    Z[1:k+1] = [r[n]^(ℓ+1) * (real(o[n]) + im * (imag(o[n]) + real(o[n]) * num/r[n])) for n=1:k+1]

    Z[1] = ℓ > 0 ? Complex{T}(0) : T(0)  + im * p

    def.pos.Na = k+1

    return Z

end
# ..............................................................................
function OUTSCH_WKB(E::T, grid::Grid{T}, def::Def{T}) where T<:Real

    N = grid.N
    r = grid.r
    k = grid.k
    v = def.pot
    s = def.scr
    n = def.pos.Nlctp

    n > 0 || error("Error: OUTSCH_WKB requires non-zero lower classical turning point")

    p = sqrt.(abs.(v .+ s .- E))                             # quasi-classical momentum
    I = [grid_trapezoidal_integral(p, i:n, grid) for i=1:n]  # quasi-classical integral
    P = exp.(-I) ./ sqrt.(p[1:n])                            # WKB solution
    P = append!(P,zeros(N-n))
    Q = grid_lagrange_derivative(P, grid; k)

    Na = def.pos.Na = findfirst(x -> abs(x) > 1.0e-10, P)

    Na > k+1 || error("Error: Na ≤ k+1 (quasi-classical approximation marginal)")

    return P .+ im * Q

end

# =========================== Adams sector =====================================

struct Adams{T}

    G::Vector{Matrix{T}}
    σ::Vector{Matrix{T}}
    Minv::Vector{Matrix{T}}
    Z::Vector{Complex{T}}

end

function castAdams(E::T, grid::Grid{T}, def::Def{T}) where T<:Real

    am = def.am
     k = def.k
     ℓ = def.orbit.ℓ

    N = def.pos.N
    def.pos.Nmin = _get_Nmin(def)
    def.pos.Nlctp = _get_Nlctp(E, def)
    def.pos.Nuctp = _get_Nuctp(E, def)

    G = matG(E, grid, def)
    σ = matσ(E, grid, def)
    M = matMinv(E, grid, def, am[end])
    Z = ℓ < 5 ? OUTSCH(grid, def, σ) : OUTSCH_WKB(E, grid, def)

    return Adams(G, σ, M, Z)

end

function updateAdams!(adams::Adams{T}, E, grid::Grid{T}, def::Def{T}) where T<:Real

    E = myconvert(def.T, E)
    G = matG(E, grid, def)
    σ = matσ(E, grid, def)
    M = matMinv(E, grid, def, def.am[end])
    Z = adams.Z

    def.pos.Na = _get_Na(Z, def)
    def.pos.Nlctp = _get_Nlctp(E, def)
    def.pos.Nuctp = _get_Nuctp(E, def)

    return Adams(G, σ, M, Z)

end

# ====================== INSCH sector ==========================================

function _nexta!(a::Vector{T}, aEnd::T, ℓ::T, σ::T, λ::T) where T <: Real

    k = T(length(a))

    one = T(1)
    two = T(2)
    add = (ℓ*(ℓ+one)-(σ-k)*(σ-k+one)) * aEnd / (two*λ*k)

    push!(a, add)

    return a

end
# ..............................................................................
function _nextb!(b::Vector{T}, aEnd::T, ℓ::T, σ::T) where T <: Real

    k = T(length(b))

    one = T(1)
    two = T(2)
    add = ((σ+k)*(σ-k+one)-ℓ*(ℓ+one)) * aEnd / (two*k)

    push!(b, add)

    return b

end
# ..............................................................................
function _nexto!(o::Vector{Complex{T}}, ρ::T, ρN::T, a::Vector{T}, b::Vector{T}, σ::T, λ::T) where T <: Real

    n = length(a)

    one = T(1)
    two = T(2)
    den = one

    c::Vector{T} = [1]

    for i = 2:n
        den *= ρ
        push!(c, one/den)
    end

    a′ = a ⋅ c
    b′ = b ⋅ c

    push!(o, (ρ/ρN)^σ * exp(-λ*(ρ-ρN)) * (a′ + b′*im))   #  α′=sum([α[s]ρ^(1-s) for s ∈ eachindex(α)], with α ∈ {a,b}

    return o

end
# ..............................................................................
function INSCH(E::T, grid::Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real

    N = def.pos.N
    r = grid.r
    k = def.k
    Z = adams.Z

    B = BigFloat

    Zc = myconvert(B, def.atom.Zc)
     ℓ = myconvert(B, def.orbit.ℓ)
     E = myconvert(B, E)

    λ = sqrt(-2E)
    σ = Zc/λ

    a::Vector{B} = [1]
    b::Vector{B} = [-λ]
    o = Complex{B}[]

    smax = ceil(Int,max(2σ,k))
    for s=1:smax
        _nextb!(b,a[end],ℓ,σ)
        _nexta!(a,a[end],ℓ,σ,λ)
    end

    ρN = myconvert(B, r[N])
    for n=0:k
        ρ = myconvert(B, r[N-k+n])
        _nexto!(o, ρ, ρN, a, b, σ, λ)  # Johnson (2.77) and (2,78)
    end

    o /= real(o[end])

    Z[N-k:N] = convert.(Complex{T}, o)

    def.pos.Nb = N - k

    return Z[N-k:N]

end

# ======================= adams_moulton_inward section =========================

function _prepend!(Z2, n, m, G, am, k, N)

    P = am[1:k] ⋅ [G[n+k-j][1,2] * imag(Z2[j+1]) for j=0:k-1]
    Q = am[1:k] ⋅ [G[n+k-j][2,1] * real(Z2[j+1]) for j=0:k-1]
    z = Z2[1] - (P + im*Q)
    z = m[n] * [real(z), -imag(z)]

    return prepend!(Z2, z[1] - im*z[2])  # change sign of derivative to negative

end
# ..............................................................................
function adams_moulton_inward(E::T, grid::Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real

    N = grid.N
    k = grid.k
    Z = adams.Z
    Nuctp = def.pos.Nuctp
    sgn = isodd(def.orbit.n′) ? -1.0 : 1.0

    Z2 = INSCH(E, grid, def, adams)

    for n=N-k-1:-1:Nuctp
        _prepend!(Z2, n, adams.Minv, adams.G, def.am, k, N)
    end

    ΔQ = imag(Z[Nuctp]) - imag(Z2[1])/ real(Z2)[1]

    Z[Nuctp:N] = sgn * Z2[1:N-Nuctp+1] / real(Z2)[1]  # set amplitude at c.t.p. to +1 or -1

    def.pos.Na = _get_Na(Z, def)
    def.pos.Nb = _get_Nb(Z, def)

    return ΔQ, Z

end

# ======================= adams_moulton_outward ================================

function adams_moulton_outward(def::Def{T}, adams::Adams{T}) where T<:Real

    am = def.am
    k = def.k

    Minv = adams.Minv
    G = adams.G
    Z = adams.Z

    N  = def.pos.N
    Na = def.pos.Na
    Nb = def.pos.Nb
    Nuctp = def.pos.Nuctp

    for n=Na:Nuctp-1
        P = am[1:k] ⋅ [G[n+1-k+j][1,2] * imag(Z[n+1-k+j]) for j=0:k-1]
        Q = am[1:k] ⋅ [G[n+1-k+j][2,1] * real(Z[n+1-k+j]) for j=0:k-1]
        z = Z[n] + (P + Q*im)
        z = Minv[n+1] * [real(z), imag(z)]
        Z[n+1] = z[1] + z[2]*im
    end

    Z1 = abs(real(Z[Nuctp]))

    Z[1:Nuctp] /= Z1         # set amplitude at c.t.p. to +1/-1 (nodes even/odd)

    def.pos.Na = _get_Na(Z, def)

    return Z

end

# ======================= adams_moulton_normalized =============================

function adams_moulton_normalized(Z::Vector{Complex{T}}, ΔQ::T, grid::Grid{T}, def::Def{T}) where T<:Real

    Na = def.pos.Na
    Nuctp = def.pos.Nuctp
    Nb = def.pos.Nb
    k = def.k

    norm = grid_trapezoidal_integral(real(Z) .^2, Na-k, Nb+k, grid)

    ΔE = ΔQ * abs(real(Z[Nuctp])) / T(2)

    return ΔE/norm, Z/sqrt(norm)

end

# ======================= count_nodes(Z, def) ==================================

function count_nodes(Z::Vector{Complex{T}}, def::Def{T}) where T<:Real

    Na = def.pos.Na
    Nb = def.pos.Nb
    c = [sign(real(Z[n])*real(Z[n+2])) for n=Na:2:Nb-2]
    c = findall(x -> x == -1, c)

    nodes = length(c)

    return nodes

end

# =============== solve_adams_moulton(E, grid, def, adams) =====================

function solve_adams_moulton(E::T, grid::Grid{T}, def::Def{T}, adams::Adams) where T<:Real

     adams = updateAdams!(adams, E, grid, def)
         Z = adams_moulton_outward(def, adams)
     ΔQ, Z = adams_moulton_inward(E, grid, def, adams)
     ΔE, Z = adams_moulton_normalized(Z, ΔQ, grid, def)

    return adams, ΔE, Z

end
