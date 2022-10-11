# ======================== OUTSCH(E, grid, def, σ) =============================

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
# ..............................................................................
function _vecOUTSCH(k::Int, matLD::Matrix{T}, p::T, q::T) where T<:Real

        v = [matLD[2,1]]

    for i=2:k  push!(v, matLD[1+i,1] * p)  end       # Johnson (2.67)
    for i=1:k  push!(v, matLD[1+i,1] * q)  end       # Johnson (2.68)

    return -v

end
# ..................................................................................
#function _update_Z!(Z::Vector{Complex{T}}, o::Vector{Complex{T}}, Na::Int, k::Int) where T<:Real
#
#    for n=Na:Na+k
#        Z[n] = o[n-Na+1]
#    end
#
#    return Z
#
#end
# ..............................................................................

@doc raw"""
    OUTSCH(E::T, grid::Grid{T}, def::Def{T}, σ::Vector{Matrix{T}}})) where T<:Real

Solution of the Schrödinger for the first ``k`` points on the `grid`, where
``k`` is the Adams-Moulton order. The WKB solution for energy `E` is used
when the WKB approximation is valid (for nonzero angular momentum at distances
below the inner classical turning point - ictp)

NB. The plot functions require CairoMakie (not included in CamiXon) to be
installed. For the code see `plotfunctions.jl` in the CamiXon.depot directory.
#### Example:
NB. `plot_wavefunction` (see `plot_functions.jl` in `CamiXon.depot`) uses
`CairoMakie`, which is not included in the `CamiXon` package.
```
Ecal, grid, def, adams = demo_hydrogen(n=1, ℓ=0)
Z = OUTSCH(Ecal, grid, def, adams.σ)
println("\nZ: standard Ansatz for wavefunction (n < Na=$(def.pos.Na)))")
    Orbital: 1s
        principal quantum number: n = 1
        radial quantum number: n′ = 0 (number of nodes in radial wavefunction)
        orbital angular momentum of valence electron: ℓ = 0
    Grid created: exponential, Float64, Rmax = 63.0 a.u., Ntot = 100, h = 0.1, r0 = 0.00286033
    Def created for hydrogen 1s on exponential grid

    Z: standard Ansatz for wavefunction (n < Na=8))

Ecal, grid, def, adams = demo_hydrogen(n=10, ℓ=5)
Z = OUTSCH(Ecal, grid, def, adams.σ);
println("\nZ: WKB Ansatz for wavefunction (n < Na=$(def.pos.Na)))")
    Orbital: 10h
        principal quantum number: n = 10
        radial quantum number: n′ = 4 (number of nodes in radial wavefunction)
        orbital angular momentum of valence electron: ℓ = 5
    Grid created: exponential, Float64, Rmax = 360.0 a.u., Ntot = 550, h = 0.0181818, r0 = 0.0163447
    Def created for hydrogen 10h on exponential grid

    Z: WKB Ansatz for wavefunction (n < Na=70))

plot_wavefunction(Z, 1:def.pos.Na, E, grid, def; reduced=true)
```
![Image](./assets/OUTSCH_H1_10h.png)
"""
function OUTSCH(E::T, grid::Grid{T}, def::Def{T}, σ::Vector{Matrix{T}}) where T<:Real

    N = grid.N
    r = grid.r
    k = grid.k
    v = def.pot
    s = def.scr
    Nlctp = def.pos.Nlctp

    o = Nlctp > 0 ? OUTSCH_WKB(grid, def, σ) : OUTSCH_WJ(grid, def, σ)

    return o

end
# ..............................................................................
function OUTSCH_hydrogenic(grid::Grid{T}, def::Def{T}, σ::Vector{Matrix{T}}) where T<:Real

        r = grid.r
        k = def.k
        N = def.pos.N
        ℓ = def.orbit.ℓ
        n′= def.orbit.n′
     Zval = def.atom.Z

    Z = zeros(Complex{T},N)

    P = real(Z)
    Q = imag(Z)

    P[1:k+1] = [r[n]^(ℓ+1) * exp(-Zval/(n′+ℓ+1)*r[n]) for n=1:k+1]
    Q[1:k+1] = [(ℓ+1)*r[n]^ℓ - r[n]^(ℓ+1)*exp(-Zval/(n′+ℓ+1)*r[n])*Zval/(n′+ℓ+1) for n=1:k+1]

    def.pos.Na = k+1

    return P .+ im * Q

end
# ..............................................................................
function OUTSCH_WKB(E::T, grid::Grid{T}, def::Def{T}, σ::Vector{Matrix{T}}) where T<:Real

    N = grid.N
    r = grid.r
    k = grid.k
    v = def.pot
    s = def.scr
    Nlctp = def.pos.Nlctp

    p = sqrt.(abs.(v .+ s .- E))                                     # quasi-classical momentum
    I = [grid_trapezoidal_integral(p, i:Nlctp, grid) for i=1:Nlctp]  # quasi-classical integral
    P = exp.(-I) ./ sqrt.(p[1:Nlctp])                                # WKB solution
    P = append!(P,zeros(N-Nlctp))
    Q = grid_differentiation(P, grid)

    Na = def.pos.Na = findfirst(x -> abs(x) > 1.0e-10, P)

    Na > k+1 || return OUTSCH_WalterJohnson(grid, def, σ)

    return P .+ im * Q

end
function OUTSCH_WJ(grid::Grid{T}, def::Def{T}, σ::Vector{Matrix{T}}) where T<:Real

        r = grid.r
        k = def.k
        N = def.pos.N
        ℓ = def.orbit.ℓ
     Zval = def.atom.Z
    matLD = def.matLD

    p = T(1)
    q = convert(T, -Zval//(ℓ + 1) )
    num = T(ℓ + 1)

    o = zeros(Complex{T}, k+1)
    o[1] = p + q * im                       # boundary condition

    m = _matOUTSCH(k, matLD, σ)
    v = _vecOUTSCH(k, matLD, p, q)

    u = inv(m) * v                          # solve 2dx2d set of equations

    for n=1:k
        o[n+1] = u[n] + im * u[n+k]
    end

    Z = zeros(Complex{T},N)

    Z[1:k+1] = [r[n]^(ℓ+1) * (real(o[n]) + im * (imag(o[n]) + real(o[n]) * num/r[n])) for n=1:k+1]

    Z[1] = ℓ > 0 ? Complex{T}(0) : T(0)  + im * p

    def.pos.Na = k+1

    return Z

end
