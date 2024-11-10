# SPDX-License-Identifier: MIT

# author: Jook Walraven - 14-2-2023

# ==============================================================================
#                               outsch.jl
# ==============================================================================

function kanweg_matOUTSCH(k::Int, matLD::Matrix{T}, matσ::Vector{Matrix{T}}) where T<:Real

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
function _matOUTSCH(k::Int, matLD::Matrix{T}, σ::Vector{Matrix{T}}) where T<:Real

    mat = Array{T,2}(undef,2k,2k)
    nul = T(0)

    for i=1:k 
        for j=1:k
            mat[i,j] = matLD[1+i,1+j] # set blok 1,1
            mat[i,k+j] = nul          # prepare blok 1,2
            mat[k+i,j] = nul          # prepare blok 2,1
            mat[k+i,k+j] = mat[i,j]   # prepare blok 2,2
        end
    end
    
    for i=1:k
        mat[i,k+i] = -σ[i][1,2]       # set blok 1,2
        mat[k+i,i] = -σ[i][2,1]       # set blok 2,1
        mat[k+i,k+i] -= σ[i][2,2]     # set blok 2,2
    end

    return mat

end
# ..............................................................................
function kanweg_vecOUTSCH(k::Int, matLD::Matrix{T}, p::T, q::T) where T<:Real

        v = [matLD[2,1]]

    for i=2:k  Base.push!(v, matLD[1+i,1] * p)  end       # Johnson (2.67)
    for i=1:k  Base.push!(v, matLD[1+i,1] * q)  end       # Johnson (2.68)

    return -v

end
function _vecOUTSCH(k::Int, matLD::Matrix{T}, p::T, q::T) where T<:Real

    v = [matLD[2,1]]

    for i=2:k  Base.push!(v, matLD[1+i,1] * p)  end       # Johnson (2.67)
    for i=1:k  Base.push!(v, matLD[1+i,1] * q)  end       # Johnson (2.68)

    return -v

end

# ------------------------------------------------------------------------------
#                       INSCH!(Z, E, grid, def)
# ------------------------------------------------------------------------------

@doc raw"""
    OUTSCH!(Z::Vector{Complex{T}}, E::T, grid::Grid{T}, def::Def{T}, adams::Adams1{T}) where T<:Real

Ansatz solution for the *outward* integration of the radial wave equation for the first ``k`` points 
on the [`Grid`](@ref), where ``k`` is the Adams-Moulton order. For angular momentum `0 ≤ ℓ ≤ 5` the 
Walter Johnson Ansatz is used; for ``ℓ > 5`` the Ansatz is based on the WKB solution for energy `E` 
at distances *far below* the inner classical turning point - ictp)

#### Example:
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

plot_wavefunction(Z, 1:def.pos.Na, grid, def; reduced=true)
```
The plot is made using `CairomMakie`.
NB.: `plot_wavefunction` is not included in the `CamiXon` package.
![Image](./assets/OUTSCH_H1_10h.png)
"""
function OUTSCH(E::T, grid::Grid{T}, def::Def{T}, σ::Vector{Matrix{T}}) where T<:Real

    N = grid.N
    r = grid.r
    k = grid.k
    v = def.pot
    s = def.scr
    Nlctp = def.pos.Nlctp

    o = Nlctp > 0 ? OUTSCH_WKB(E, grid, def, σ) : OUTSCH_WJ(grid, def, σ)

    return o

end
function OUTSCH!(Z::Vector{Complex{T}}, E::T, grid::Grid{T}, def::Def{T}, adams::Adams1{T}) where T<:Real

    ℓ = def.orbit.ℓ
    
    Z = ℓ < 6 ? OUTSCH_WJ!(Z, grid, def, adams) : OUTSCH_WKB!(Z, E, grid, def)

    return Z

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
function OUTSCH_WJ!(Z::Vector{Complex{T}}, grid::Grid{T}, def::Def{T}, adams::Adams1{T}) where T<:Real

    r = grid.r
    k = def.k
    N = def.pos.N
    ℓ = def.orbit.ℓ
 Zval = def.atom.Z
matLD = def.matLD
    σ = adams.σ

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

    Z[1:k+1] = [r[n]^(ℓ+1) * (real(o[n]) + im * (imag(o[n]) + real(o[n]) * num/r[n])) for n=1:k+1]

    Z[1] = ℓ > 0 ? Complex{T}(0) : T(0)  + im * p

    def.pos.Na = k+1

    return Z

end

# ..............................................................................
function OUTSCH_WKB(E::T, grid::Grid{T}, def::Def{T}, σ::Vector{Matrix{T}}) where T<:Real

    N = grid.N
    r = grid.r
    k = grid.k
    v = def.pot
    s = def.scr
    cWKB = def.pos.cWKB
    Nlctp = def.pos.Nlctp

    p = sqrt.(abs.(v .+ s .- E))                                     # quasi-classical momentum
    I = [grid_integration(p, i:Nlctp, grid) for i=1:Nlctp]           # quasi-classical integral
    P = exp.(-I) ./ sqrt.(p[1:Nlctp])                                # WKB solution
    P = append!(P,zeros(N-Nlctp))
    Q = grid_differentiation(P, grid)

    Na = def.pos.Na = findfirst(x -> abs(x) > cWKB, P)

    Na > k+1 || return OUTSCH_WJ(grid, def, σ)

    return P .+ im * Q

end
function OUTSCH_WKB!(Z::Vector{Complex{T}}, E::T, grid::Grid{T}, def::Def{T}) where T<:Real
    
    N = grid.N
    ################ r = grid.r
    ################ r′= grid.r′
    k = grid.k
    Nlctp = def.pos.Nlctp
    ################ two = T(2)
    pot = def.potscr

    p = Array{T,1}(undef,N)
    I = Array{T,1}(undef,N)
    P = Array{T,1}(undef,N)
    Q = Array{T,1}(undef,N)

    for n=Nlctp:-1:1
        p[n] = sqrt(abs(pot[n]-E))                             # quasi-classical momentum
        I[n] = grid_integration(p, n:Nlctp, grid)
        P[n] = exp(-I[n])/sqrt(p[n])                           # WKB solution
        def.pos.Na = P[n] > 1.0e-30 ? n : break
    end

    Na = def.pos.Na = def.pos.Na + k
    
    Q[Na-k:Na+k] = grid_differentiation1(P, Na-k:Na+k, grid)   # avoid lower end point correction by doubling range
    
    for n=Na-k:Na+k
        Z[n] = P[n] + im * Q[n]
    end
    
    return Z

end