
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

