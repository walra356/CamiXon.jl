# ======================== Pos(Nmin, Na, Nuctp, Nb, N, nodes) ===================


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

    N = grid.N
    r = grid.r
    k = grid.k
    v = def.pot
    s = def.scr
    cWKB = def.pos.cWKB
    Nuctp = def.pos.Nuctp

    p = sqrt.(abs.(v .+ s .- E))                            # quasi-classical momentum
    I = [grid_integration(p, Nuctp:i, grid) for i=Nuctp:N]  # quasi-classical integral
    P = exp.(-I) ./ sqrt.(p[Nuctp:N])                       # WKB solution
    P = prepend!(P,zeros(Nuctp-1))
    Q = grid_differentiation(P, grid)

    Nb = findfirst(x -> 0 < abs(x) < cWKB, P)
    Nb = def.pos.Nb = isnothing(Nb) ? N-k : min(N-k, Nb)

    Z = P .+ im * Q

    Z2 = copy(Z[Nb:N])

    return Z2

end
function INSCH_WJ(E::T, grid::Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real

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
