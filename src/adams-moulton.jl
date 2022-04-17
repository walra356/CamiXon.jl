"""
    matG(E::T, grid::Grid{T}, def::Def{T}) where T<:Real

"""
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
@doc raw"""
    matσ(E::T, grid::Grid{T}, def::Def{T}) where T<:Real

"""
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
@doc raw"""
    matMinv(E::T, grid::Grid{T}, def::Def{T}, amEnd::T) where T<:Real

"""
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
# ..............................................................................

@doc raw"""
    OUTSCH(grid::Grid{T}, def::Def{T}, σ::Vector{Matrix{T}}) where T<:Real

"""
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

@doc raw"""
    OUTSCH_WKB(E::T, grid::Grid{T}, def::Def{T}) where T<:Real

"""
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

@doc raw"""
    Adams

* G: (`:Vector{Matrix{T}}`)
* σ: (`:Vector{Matrix{T}}`)
* Minv: (`:Vector{Matrix{T}}`)
* Z: (`:Vector{Complex{T}}`)
"""
struct Adams{T}

    G::Vector{Matrix{T}}
    σ::Vector{Matrix{T}}
    Minv::Vector{Matrix{T}}
    Z::Vector{Complex{T}}

end

@doc raw"""
    castAdams(E::T, grid::Grid{T}, def::Def{T}) where T<:Real

"""
function castAdams(E::T, grid::Grid{T}, def::Def{T}) where T<:Real

    am = def.am
     k = def.k
     ℓ = def.orbit.ℓ

    N = def.pos.N
    def.pos.Nmin = get_Nmin(def)
    def.pos.Nlctp = get_Nlctp(E, def)
    def.pos.Nuctp = get_Nuctp(E, def)

    G = matG(E, grid, def)
    σ = matσ(E, grid, def)
    M = matMinv(E, grid, def, am[end])
    Z = ℓ < 5 ? OUTSCH(grid, def, σ) : OUTSCH_WKB(E, grid, def)

    return Adams(G, σ, M, Z)

end

@doc raw"""
    updateAdams!(adams::Adams{T}, E, grid::Grid{T}, def::Def{T}) where T<:Real

"""
function updateAdams!(adams::Adams{T}, E, grid::Grid{T}, def::Def{T}) where T<:Real

    E = myconvert(def.T, E)
    G = matG(E, grid, def)
    σ = matσ(E, grid, def)
    M = matMinv(E, grid, def, def.am[end])
    Z = adams.Z

    def.pos.Na = get_Na(Z, def)
    def.pos.Nlctp = get_Nlctp(E, def)
    def.pos.Nuctp = get_Nuctp(E, def)

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

@doc raw"""
    INSCH(E::T, grid::Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real

"""
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

@doc raw"""
    adams_moulton_inward(E::T, grid::Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real

"""
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

    def.pos.Na = get_Na(Z, def)
    def.pos.Nb = get_Nb(Z, def)

    return ΔQ, Z

end

# ======================= adams_moulton_outward ================================

@doc raw"""
    adams_moulton_outward(def::Def{T}, adams::Adams{T}) where T<:Real

"""
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

    def.pos.Na = get_Na(Z, def)

    return Z

end

# ======================= adams_moulton_normalized =============================

@doc raw"""
    adams_moulton_normalized(Z::Vector{Complex{T}}, ΔQ::T, grid::Grid{T}, def::Def{T}) where T<:Real

"""
function adams_moulton_normalized(Z::Vector{Complex{T}}, ΔQ::T, grid::Grid{T}, def::Def{T}) where T<:Real

    Na = def.pos.Na
    Nuctp = def.pos.Nuctp
    Nb = def.pos.Nb
    k = def.k

    norm = grid_trapezoidal_integral(real(Z) .^2, Na-k, Nb+k, grid)

    ΔE = ΔQ * abs(real(Z[Nuctp])) / T(2)

    return ΔE/norm, Z/sqrt(norm)

end

# =============== solve_adams_moulton(E, grid, def, adams) =====================

@doc raw"""
    solve_adams_moulton(E::T, grid::Grid{T}, def::Def{T}, adams::Adams) where T<:Real

"""
function solve_adams_moulton(E::T, grid::Grid{T}, def::Def{T}, adams::Adams) where T<:Real

     adams = updateAdams!(adams, E, grid, def)
         Z = adams_moulton_outward(def, adams)
     ΔQ, Z = adams_moulton_inward(E, grid, def, adams)
     ΔE, Z = adams_moulton_normalized(Z, ΔQ, grid, def)

    return adams, ΔE, Z

end
