# SPDX-License-Identifier: MIT

# author: Jook Walraven - 14-2-2023

# ==============================================================================
#                               insch.jl
# ==============================================================================

# ======================== Pos(Nmin, Na, Nuctp, Nb, N, nodes) ===================


# ====================== INSCH sector ==========================================

function _nexta!(a::Vector{T}, aEnd::T, ℓ::T, σ::T, λ::T) where T <: Real

    k = T(length(a))

    one = T(1)
    two = T(2)
    add = (ℓ*(ℓ+one)-(σ-k)*(σ-k+one)) * aEnd / (two*λ*k)

    Base.push!(a, add)

    return a

end
# ..............................................................................
function _nextb!(b::Vector{T}, aEnd::T, ℓ::T, σ::T) where T <: Real

    k = T(length(b))

    one = T(1)
    two = T(2)
    add = ((σ+k)*(σ-k+one)-ℓ*(ℓ+one)) * aEnd / (two*k)

    Base.push!(b, add)

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
        Base.push!(c, one/den)
    end

    a′ = a ⋅ c
    b′ = b ⋅ c

    Base.push!(o, (ρ/ρN)^σ * exp(-λ*(ρ-ρN)) * (a′ + b′*im))   #  α′=sum([α[s]ρ^(1-s) for s ∈ eachindex(α)], with α ∈ {a,b}

    return o

end

# ------------------------------------------------------------------------------
#                       INSCH!(Z, E, grid, def)
# ------------------------------------------------------------------------------

@doc raw"""
    INSCH(E::T, grid::Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real
 
Ansatz solution for the *inward* integration of the radial wave equation for the first ``k`` points 
on the [`Grid`](@ref), where ``k`` is the Adams-Moulton order. The Ansatz is based on the WKB solution 
for energy `E` at distances *far above* the upper classical turning point - uctp)
"""
function INSCH(E::T, grid::Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real

    N = grid.N
    r = grid.r
    k = grid.k
    v = def.pot
    s = def.scr
    Nuctp = def.pos.Nuctp

    p = sqrt.(abs.(v .+ s .- E))                            # quasi-classical momentum
    I = [grid_integration(p, Nuctp:i, grid) for i=Nuctp:N]  # quasi-classical integral
    P = exp.(-I) ./ sqrt.(p[Nuctp:N])                       # WKB solution
    P = prepend!(P,zeros(Nuctp-1))
    Q = grid_differentiation(P, grid)

    #Nb = findfirst(x -> 0 < abs(x) < cWKB, P)

    Nb = findlast(x -> abs(x) > 1.0e-7, P)
    Nb = def.pos.Nb = isnothing(Nb) ? N-k : min(N-k, Nb)

    Z = P .+ im * Q

    Z2 = copy(Z[Nb:N])

    return Z2

end

@doc raw"""
    INSCH!(Z::Vector{Complex{T}}, E::T, grid::Grid{T}, def::Def{T}) where T<:Real
 
Ansatz solution for the *inward* integration of the radial wave equation for the first ``k`` points 
on the [`Grid`](@ref), where ``k`` is the Adams-Moulton order. The Ansatz is based on the WKB solution 
for energy `E` at distances *far above* the upper classical turning point - uctp)
"""
function INSCH!(Z::Vector{Complex{T}}, E::T, grid::Grid{T}, def::Def{T}) where T<:Real
    
    Z = INSCH_WKB!(Z, E, grid, def)

    return Z

end

# ------------------------------------------------------------------------------
#                       INSCH_WKB!(Z, E, grid, def)
# ------------------------------------------------------------------------------

@doc raw"""
    INSCH_WKB!(Z::Vector{Complex{T}}, E::T, grid::Grid{T}, def::Def{T}) where T<:Real

WKB Ansatz of `k+1` points for [`INSCH!`](@ref)
"""
function INSCH_WKB!(Z::Vector{Complex{T}}, E::T, grid::Grid{T}, def::Def{T}) where T<:Real

    N = grid.N
    r = grid.r
    k = grid.k
    Nuctp = def.pos.Nuctp
    two = T(2)
    pot = def.potscr

    p = Array{T,1}(undef,N)
    I = Array{T,1}(undef,N)
    P = Array{T,1}(undef,N)
    Q = Array{T,1}(undef,N)
    
    for n=Nuctp:N
        p[n] = sqrt(abs(pot[n]-E))                         # quasi-classical momentum   
        I[n] = grid_integration(p, Nuctp:n, grid)
        P[n] = exp(-I[n])/sqrt(p[n])                       # WKB solution
        def.pos.Nb = P[n] > 1.0e-30 ? n : break
    end

    Nb = def.pos.Nb = def.pos.Nb - k
    
    Q[Nb-k:Nb+k] = grid_differentiation(P, grid, Nb-k:Nb+k)   # avoid lower end point correction by doubling range
    
    for n=Nb-k:Nb+k
        Z[n] = P[n] + im * Q[n]
    end 
    
    return Z

end

# ------------------------------------------------------------------------------
#                       INSCH_WJ!(Z, E, grid, def)
# ------------------------------------------------------------------------------

function INSCH_WJ(E::T, grid::Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real

############ "kept for the record" #####################

    N = def.pos.N
    r = grid.r
    k = def.k
    Z = adams.Z

    B = BigFloat

    Zc = convert(B, def.atom.Zc)
     ℓ = convert(B, def.orbit.ℓ)
     E = convert(B, E)

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

    ρN = convert(B, r[N])
    for n=0:k
        ρ = convert(B, r[N-k+n])
        _nexto!(o, ρ, ρN, a, b, σ, λ)  # Johnson (2.77) and (2,78)
    end

    o /= real(o[end])

    Z[N-k:N] = convert.(Complex{T}, o)

    def.pos.Nb = N - k

    return Z[N-k:N]

end
function INSCH_WJ!(Z::Vector{Complex{T}}, E::T, grid::Grid{T}, def::Def{T}, adams::Adams1{T}, a::Vector{BigFloat}, b::Vector{BigFloat}, c::Vector{BigFloat}, o::Vector{Complex{BigFloat}}) where T<:Real

    ############ "kept for the record" #####################
    
    N = grid.N
    r = grid.r
    k = def.k
    
    B = BigFloat

    Zc = convert(B, def.atom.Zc)
     ℓ = convert(B, def.orbit.ℓ)
     E = convert(B, E)

    one = T(1)
    two = T(2)
    den = one

    λ = sqrt(-2E)
    ζ = Zc/λ
    
    a[1] = T(1)
    b[1] = -λ
    c[1] = T(1)

    smax = ceil(Int,max(2ζ,k))
    for s=2:smax
        b[s] = ((ζ+k)*(ζ-k+one)-ℓ*(ℓ+one)) * a[s-1] / (two*k)
        a[s] = (ℓ*(ℓ+one)-(ζ-k)*(ζ-k+one)) * a[s-1] / (two*λ*k)
    end

    a1 = a[1:smax]
    b1 = b[1:smax]
    c1 = c[1:smax]
    
    ρN = convert(B, r[N])
    for n=0:k
        den = one
        ρ = convert(B, r[N-k+n])
        for i=2:smax
            den *= ρ
            c1[i] = one/den
        end
        a′ = a1 ⋅ c1
        b′ = b1 ⋅ c1
        o[n+1] = (ρ/ρN)^ζ * exp(-λ*(ρ-ρN)) * (a′ + b′*im) # Johnson (2.77) and (2,78)
    end

    #o .*= 1.0e-300 #real(o[1])

    Z[N-k:N] = convert.(Complex{T}, o)

    def.pos.Nb = N - k

    return Z

end

