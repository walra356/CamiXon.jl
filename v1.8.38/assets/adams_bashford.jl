# SPDX-License-Identifier: MIT

# author: Jook Walraven - 4-9-2024

# ==============================================================================
#                         adams-bashford.jl
# ==============================================================================


# ------------------------------------------------------------------------------
#                       adams_bashford_outward(Z, def, adams)
# ------------------------------------------------------------------------------

function adams_bashford_outward(Z::Vector{Complex{T}}, def::Def{T}, adams::Adams{T}) where T<:Real

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

    norm = abs(real(Z[Nuctp]))

    Z[1:Nuctp] /= norm    # set amplitude at u.c.t.p. to +1/-1 (nodes even/odd)

    #def.pos.Na = get_Na(Z, def)

    return Z

end