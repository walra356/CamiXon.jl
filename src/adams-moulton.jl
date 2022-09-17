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

# ======================= adams_moulton_inward section =========================

function _prepend!(Z2, n, m, G, am, k)

    P = am[1:k] ⋅ [G[n+k-j][1,2] * imag(Z2[k-j]) for j=0:k-1]
    Q = am[1:k] ⋅ [G[n+k-j][2,1] * real(Z2[k-j]) for j=0:k-1]
    z = Z2[1] - (P + im*Q)
    z = m[n] * [real(z), -imag(z)]

    return prepend!(Z2, z[1] - im*z[2])               # change sign of derivative to negative

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
        _prepend!(Z2, n, adams.Minv, adams.G, def.am, k)
    end

    ΔQ = imag(Z[Nuctp]) - imag(Z2[1])/ real(Z2)[1]

    Z[Nuctp:N] = sgn * Z2[1:N-Nuctp+1] / real(Z2)[1]  # set amplitude at c.t.p. to +1 or -1

    def.pos.Na = get_Na(Z, def)
    def.pos.Nb = get_Nb(Z, def)

    return ΔQ, Z

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

# =============== adams_moulton_solve(E, grid, def, adams) =====================

@doc raw"""
    adams_moulton_solve(E::T, grid::Grid{T}, def::Def{T}, adams::Adams) where T<:Real

"""
function adams_moulton_solve(E::T, grid::Grid{T}, def::Def{T}, adams::Adams) where T<:Real

     adams = updateAdams!(adams, E, grid, def)
         Z = adams_moulton_outward(def, adams)
     ΔQ, Z = adams_moulton_inward(E, grid, def, adams)
     ΔE, Z = adams_moulton_normalized(Z, ΔQ, grid, def)

    return adams, ΔE, Z

end
