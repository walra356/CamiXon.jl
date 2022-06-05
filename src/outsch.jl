# ======================== OUTSCH(grid, def, σ) ================================


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
    q = convert(T, -Zval//(ℓ + 1) )
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
    Q = grid_differentiation(P, grid; k)

    Na = def.pos.Na = findfirst(x -> abs(x) > 1.0e-10, P)

    Na > k+1 || error("Error: Na ≤ k+1 (quasi-classical approximation marginal)")

    return P .+ im * Q

end
