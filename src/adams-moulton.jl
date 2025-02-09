# SPDX-License-Identifier: MIT

# author: Jook Walraven - 4-9-2024

# ==============================================================================
#                         adams-moulton.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#                       matG(E, grid, def) 
# ------------------------------------------------------------------------------

"""
    matG(E::T, pot::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real

coupling matrix - Johnson (2.54)
"""
function matG(E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real
# ==============================================================================
# matG - coupling matrix - Johnson (2.54)
# ==============================================================================
    
    G = def.G
    r′= grid.r′

    two = T(2)
    nul = T(0)
    pot = def.potscr

    for n ∈ eachindex(G)
        G[n][1,1] = nul
        G[n][1,2] =  r′[n]
        G[n][2,1] =  r′[n] * two * (-E + pot[n])
        G[n][2,2] = nul
    end

    return G

end

# ------------------------------------------------------------------------------
#                       matσ(E, grid, def) 
# ------------------------------------------------------------------------------

@doc raw"""
    matσ(E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real

coupling matrix - Johnson (2.54)
"""
function matσ(E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real
# ==============================================================================
# matσ - coupling matrix - Johnson (2.54)
# ==============================================================================
    
    σ = def.σ
    
    r = grid.r
    r′= grid.r′
    Zval = def.atom.Z
    ℓ = def.orbit.ℓ
    s = def.scr

    nul = T(0)
    one = T(1)
    two = T(2)
    num = -T(ℓ + 1)

    for n ∈ eachindex(σ)
        σ[n][1,1] = nul
        σ[n][1,2] = r′[n]                                  # b in Johnson (2.69)
        σ[n][2,1] = r′[n] * two * (-E - Zval/r[n] + s[n])  # c in Johnson (2.69)
        σ[n][2,2] = r′[n] * two * num / r[n]               # d in Johnson (2.69)
    end

    r1 = T(eps(Float64))  # quasi zero
    σ[1][2,1] = r′[1] * two * (-E - Zval / r1 + s[1])
    σ[1][2,2] = r′[1] * two * num / r1

    return σ

end

# ------------------------------------------------------------------------------
#                       matMinv(E, grid, def, amEnd) 
# ------------------------------------------------------------------------------

@doc raw"""
    matMinv(E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real

Adams-Moulton correction matrix - Johnson (2.56)
"""
function matMinv(E::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real
# ==============================================================================
# matMinv - Adams-Moulton correction matrix - Johnson (2.56)
# ==============================================================================

    Minv = def.Minv
    r′= grid.r′
    amEnd = def.am[end]
    
    one = T(1)
    two = T(2)
    pot = def.potscr
    
    for n ∈ eachindex(Minv) # Johnson (2.56)
        Minv[n][1,1] = one
        Minv[n][1,2] = r′[n] * amEnd
        Minv[n][2,1] = r′[n] * amEnd * two * (-E + pot[n])
        Minv[n][2,2] = one
    end
    
    for n ∈ eachindex(Minv)
        Minv[n] ./= (one - Minv[n][1,2]*Minv[n][2,1])
    end

    return Minv

end

# ------------------------------------------------------------------------------
#                  adams_moulton_outward!(def, adams)
# ------------------------------------------------------------------------------

function _updatenodes!(nodes::Int, z::T) where T<:Real

    s = iseven(nodes) ? true : false
    
    if (z < T(0)) == s nodes += 1 end

    return nodes
        
end

@doc raw"""
    adams_moulton_outward!(def::Def{T}, adams::Adams{T}) where T<:Real

"""
function adams_moulton_outward!(Z::Vector{Complex{T}}, def::Def{T}, adams::Adams{T}) where T<:Real

    k = def.k
    am = def.am[1:k]

    P = Array{T,1}(undef,k)
    Q = Array{T,1}(undef,k)

    m = adams.Minv
    G = adams.G

    Na = def.pos.Na
    Nuctp = def.pos.Nuctp
    Nlctp = def.pos.Nlctp

    nodes = 0
   
    for n=Na:Nuctp-1
        for j=0:(k-1)
            P[j+1] = G[n+1-k+j][1,2] * imag(Z[n+1-k+j])
            Q[j+1] = G[n+1-k+j][2,1] * real(Z[n+1-k+j])
        end
        z = Z[n] + (LinearAlgebra.dot(am, P) + im*LinearAlgebra.dot(am, Q))
        z1 = m[n+1][1,1] * real(z) + m[n+1][1,2] * imag(z)          
        z2 = m[n+1][2,1] * real(z) + m[n+1][2,2] * imag(z)
        Z[n+1] = z1 + im*z2
        nodes = _updatenodes!(nodes, z1)
    end

    def.pos.nodes = nodes

    num = Nlctp > 0 ? sign(real(Z[Nlctp])) : T(1.0)
    
    norm = num / abs(real(Z[Nuctp])) 

    for n=1:Nuctp
        Z[n] *= norm
    end

    return Z

end

# ------------------------------------------------------------------------------
#                  adams_moulton_inward!(E, grid, def, adams)
# ------------------------------------------------------------------------------

@doc raw"""
    adams_moulton_inward!(E::T, grid::CamiDiff.Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real

"""
function adams_moulton_inward!(Z::Vector{Complex{T}}, def::Def{T}, adams::Adams{T}) where T<:Real

    k = def.k
    am = def.am[1:k]

    amP = Array{T,1}(undef,k)
    amQ = Array{T,1}(undef,k)

    m = adams.Minv
    G = adams.G

    N  = def.pos.N
    Nb = def.pos.Nb
    Nuctp = def.pos.Nuctp
    
    P0 = real(Z[Nuctp])
    Q0 = imag(Z[Nuctp])
    
    Nb+k ≤ N || error("Error: grid boost required (increase Ntot)")

    for n=Nb:-1:Nuctp
        for j=0:(k-1)
            amP[j+1] = G[n+k-j][1,2] * imag(Z[n+k-j])
            amQ[j+1] = G[n+k-j][2,1] * real(Z[n+k-j])
        end
        z = Z[n+1] - (LinearAlgebra.dot(am, amP) + im*LinearAlgebra.dot(am, amQ))
        z1 = m[n][1,1] * real(z) - m[n][1,2] * imag(z)          
        z2 = m[n][2,1] * real(z) - m[n][2,2] * imag(z)
        Z[n] = z1 - im*z2
    end
   
    norm = P0/real(Z[Nuctp])
    
    for n=Nuctp:N
        Z[n] *= norm
    end
    
    ΔQ = imag(Z[Nuctp]) - Q0

    return ΔQ, Z

end

# ------------------------------------------------------------------------------
#                  adams_moulton_normalize!(Z, ΔQ, grid, def)
# ------------------------------------------------------------------------------

@doc raw"""
    adams_moulton_normalize!(Z::Vector{Complex{T}}, ΔQ::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real

"""
function adams_moulton_normalize!(Z::Vector{Complex{T}}, ΔQ::T, grid::CamiDiff.Grid{T}, def::Def{T}) where T<:Real

    Nuctp = def.pos.Nuctp
    
    N = grid.N

    norm = CamiDiff.grid_integration(real(Z) .^2, grid, 1, N)

    ΔE = ΔQ * abs(real(Z[Nuctp])) / T(2)

    ΔE = ΔE/norm
    Z = Z/sqrt(norm)

    return ΔE, Z

end

# ------------------------------------------------------------------------------
#                  adams_moulton_solve!(Z, E, grid, def, adams)
# ------------------------------------------------------------------------------

@doc raw"""
    adams_moulton_solve!(Z::Vector{Complex{T}}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}, adams::Adams1{T}) where T<:Real

Numerical solution of the 1D Schrödinger equation for the radial motion of a
*valence* electron of energy `E`. Output: the improved `Adams` object, the
energy convergence `ΔE`, and `Z`, where `P = real(Z)` is the *reduced* radial
wavefunction and `Q = imag(Z)` its derivative.
#### Example:
```
atom = castAtom(Z=1, A=1, Q=0, msg=true)
orbit = castOrbit(n=1, ℓ=0)
grid = autoGrid(atom, orbit, Float64; Nboost=1, msg=true)
def = castDef(grid, atom, orbit, codata)
E = Ecal = convert(grid.T, bohrformula(atom.Z, orbit.n))
adams = castAdams(E, grid, def);

adams, ΔE, Z = adams_moulton_solve(E, grid, def, adams)
plot_wavefunction(Z, 1:grid.N, grid, def; reduced=true)
```
The plot is made using CairomMakie.
NB.: `plot_wavefunction` is not part of the `CamiXon` package.
![Image](./assets/hydrogen-1s-prepared.png)
"""
function adams_moulton_solve!(Z::Vector{Complex{T}}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real

    updatePos!(def.pos, E, def.potscr, grid)
    adams = updateAdams!(adams, E, grid, def)
        
    Z = OUTSCH!(Z, E, grid, def, adams)
    Z = adams_moulton_outward!(Z, def, adams)
    Z = INSCH!(Z, E, grid, def)
ΔQ, Z = adams_moulton_inward!(Z, def, adams)
ΔE, Z = adams_moulton_normalize!(Z, ΔQ, grid, def)

    return adams, ΔE, Z

end
@doc raw"""
    adams_moulton_solve_refine!(Z::Vector{Complex{T}}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real
    
"""
function adams_moulton_solve_refine!(Z::Vector{Complex{T}}, E::T, grid::CamiDiff.Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real

    adams = updateAdams!(adams, E, grid, def)
    
    Z = adams_moulton_outward!(Z, def, adams)
    
ΔQ, Z = adams_moulton_inward!(Z, def, adams)

ΔE, Z = adams_moulton_normalize!(Z, ΔQ, grid, def)

    return adams, ΔE, Z

end

# -----------------------------------------------------------------------------------------
#                 adams_moulton_nodes(E, scr, grid, def; imax=25, msg=true)
# -----------------------------------------------------------------------------------------

function _strΔt(tstop::T, tstart::T) where T<:Real

    Δt = tstop-tstart
    str = Δt > 1.0  ? (repr(Δt, context=:compact => true) * " sec")      :
          Δt > 1e-3 ? (repr(Δt*1e3, context=:compact => true) * " msec") :
                      (repr(Δt*1e6, context=:compact => true) * " μsec")
    return str

end
# ..........................................................................................................
@doc raw"""
    adams_moulton_nodes(E::Real, scr::Vector{T}, grid::CamiDiff.Grid{T}, def::Def{T}; imax=25, msg=true) where T<:Real
    
"""
function adams_moulton_nodes(E::Real, scr::Vector{T}, grid::CamiDiff.Grid{T}, def::Def{T}; imax=25, msg=true) where T<:Real

    t1 = time()

    E = T(E)

    n′= def.orbit.n′     # radial quantum number (number of nodes)
    
    for n ∈ eachindex(def.pot)
        def.scr[n] = scr[n]
        def.potscr[n] = def.pot[n] + scr[n]
    end
    
    init = castInit(E, grid, def)     # init =(E[Nlctp], E[Nmin], E[Nuctp], ΔE=0.0)

    adams = castAdams(init.E, grid, def)
    
    Z = zeros(Complex{grid.T}, grid.N)  

    adams, ΔE, Z = adams_moulton_solve!(Z, init.E, grid, def, adams)
    nodes = def.pos.nodes
    msg && println("\n--- adams_moulton_nodes!:\nstart solution: $(nodes) nodes - init = ($(init.Emin), $(init.E), $(init.Emax), $(init.ΔE))")
    
    i = 0
    while def.pos.nodes ≠ n′
        adams, ΔE, Z = adams_moulton_solve!(Z, init.E, grid, def, adams) # ΔE not used for nodes ≠ n′
        init!(init, T(0), def)
        i += 1
        i < imax || break
    end
    
    nodes = def.pos.nodes

    init.ΔE = ΔE         # ΔE from here on in use (nodes == n′)
    msg && println("   initiate ΔE: $(nodes) nodes - init = ($(init.Emin), $(init.E), $(init.Emax), $(init.ΔE))")
    
    t2 = time()
    adams_moulton_report(i, init, grid, def, _strΔt(t2,t1); unitIn="Hartree", name="adams_moulton_nodes!", msg=true)
    
    nodes == n′ || error("Error: after $i iterations nodes = $(nodes) - should be $(n′) (increase imax?)")

    return def, adams, init, Z

end

# --------------------------------------------------------------------------------------------------------------
#            adams_moulton_iterate!(Z, init, grid, def, adams; imax=25, ϵ=1e-6, msg=true)
# --------------------------------------------------------------------------------------------------------------

@doc raw"""
    adams_moulton_iterate!(Z::Vector{Complex{T}}, init::Init{T}, grid::CamiDiff.Grid{T}, def::Def{T}, adams::Adams{T}; imax=25, ϵ=1e-6, msg=true) where T<:Real
    
"""
function adams_moulton_iterate!(Z::Vector{Complex{T}}, init::Init{T}, grid::CamiDiff.Grid{T}, def::Def{T}, adams::Adams{T}; imax=25, ϵ=1.0e-6, msg=true) where T<:Real
    
    t1 = time()

    n′= def.orbit.n′     # radial quantum number (number of nodes)
    nodes = def.pos.nodes
    two = T(2)
    pot = def.potscr
    
    msg && println("\n--- adams_moulton_iterate!:\nresume solution: $(nodes) nodes - init = ($(init.Emin), $(init.E), $(init.Emax), $(init.ΔE))") 
    
    Emin = init.Emin
    E = init.E
    Emax = init.Emax
    ΔE = 1.1e-3 * E

   (Emin ≤ E ≤ Emax) || error("Error: Emin ≤ E ≤ Emax interval violated")     
 
    i = 0
    while abs(ΔE/E) > 1e-3 # convergence goal
        adams, ΔE, Z = adams_moulton_solve!(Z, E, grid, def, adams)
        init = init!(init, ΔE, def)
        E = init.E
        i += 1
        i < imax || break
    end
    while abs(ΔE/E) > ϵ # convergence goal
        adams, ΔE, Z = adams_moulton_solve_refine!(Z, E, grid, def, adams)
        init = init!(init, ΔE, def)
        E = init.E
        i += 1
        i < imax || break
    end

    msg && println(                        "     upgrade ΔE: $(nodes) nodes - init = ($(init.Emin), $(init.E), $(init.Emax), $(init.ΔE))")

    str = T == Float64 ? "Float64 " : "Bigfloat "
    
    t2 = time()
    adams_moulton_report(i, init, grid, def, _strΔt(t2,t1); unitIn="Hartree", name="adams_moulton_iterate!", msg=true)
    if i == imax
        println("Warning: convergence goal not reached (increase imax and/or precision)")
    end
    
    return def, adams, init, Z

end

# --------------------------------------------------------------------------------------------------------------
#                 adams_moulton_precise!(Z, init, grid, def; imax=10, ϵ=1e-6, msg=false)
# --------------------------------------------------------------------------------------------------------------

@doc raw"""
    adams_moulton_precise!(Z, init, grid, def; imax=10, ϵ=1e-6, msg=false)
    
"""
function adams_moulton_precise!(Z, init, grid, def; imax=25, ϵ=1e-6, msg=false)

    println("\nReset parameters to BigFloat precision:")

    B = BigFloat
    Z = convert.(Complex{B}, Z)
    k = grid.k
    scr = convert.(B, def.scr)
    Rmax = convert(B, grid.r[grid.N])
    init = Init(B(init.Emin), B(init.E), B(init.Emax), B(init.ΔE))
   
    msg && println("Rmax = $(Rmax), init.E = $(init.E)")

    grid = autoGrid(def.atom, def.orbit, B; Ntot=grid.N , Rmax, k, msg)

    def = castDef(grid, def.atom, def.orbit, def.codata; def.pos, scr)
    adams = castAdams(init.E, grid, def)
    def, adams, init, Z = adams_moulton_iterate!(Z, init, grid, def, adams; imax, ϵ, msg)

    return grid, def, adams, init, Z
    
end

# --------------------------------------------------------------------------------------------------------------
#          adams_moulton_report(E, ΔE, grid, def; unitIn="Hartree", name="name" , msg=true)
# --------------------------------------------------------------------------------------------------------------

@doc raw"""
    adams_moulton_report(E::T, ΔE::T, grid::CamiDiff.Grid{T}, def::Def{T}; unitIn="Hartree", name="name" , msg=true) where T<:Real

"""
function adams_moulton_report(it::Int, init::Init{T}, grid::CamiDiff.Grid{T}, def::Def{T}, strΔT::String; unitIn="Hartree", name="name" , msg=true) where T<:Real

    ΔE = init.ΔE
    Δf = convertUnit(abs(ΔE), def.codata)
    strΔf = " (" * strValue(Δf) * ")"
    strΔE = repr(ΔE, context=:compact => true)
    strΔErel = repr(ΔE/init.E, context=:compact => true)

    str = "\nconvergence report for " * _defspecs(grid, def) * "\n"
    str *= "$T result with $(def.pos.nodes) nodes obtained after $(it) iterations in " * strΔT * "\n"
    str *= Printf.@sprintf "    binding energy: E = %.17g %s \n" init.E unitIn
    str *= ΔE ≠ 0 ? "absolute precision: ΔE = " * strΔE * " " * unitIn * strΔf * "\n" :
                    "absolute precision: ΔE = 0 (exact under $T precision)\n"
    str *= ΔE ≠ 0 ? "relative precision: ΔE/E = " * strΔErel * ""                   :
                    "relative precision: ΔE/E = 0 (exact under $T precision)"

    return msg ? println(str) : str

end

function test_adams_moulton(E::Real, scr::Vector{T}, grid::CamiDiff.Grid{T}, def::Def{T}; test=5, msg=false) where T<:Real
    
    msg && println("\n===== enter adams_moulton_nodes! =====") 
    t1 = time()

    E = T(E)

    n′= def.orbit.n′     # radial quantum number (number of nodes)
    
    for n ∈ eachindex(def.pot)
        def.scr[n] = scr[n]
        def.potscr[n] = def.pot[n] + scr[n]
    end
    
    init = castInit(E, grid, def) # init = (E[Nlctp], E[Nmin], E[Nuctp], ΔE=0.0)
    msg && println("after castInit!")  

    adams = castAdams(init.E, grid, def)
    msg && println("after castAdams")  
    
    Z = zeros(Complex{grid.T}, grid.N) 

    updatePos!(def.pos, init.E, def.potscr, grid)
    adams = updateAdams!(adams, init.E, grid, def)
    msg && println("after updateAdams")  
    msg && print("grid special points: Na = $(def.pos.Na), Nlctp = $(def.pos.Nlctp), Nmin = $(def.pos.Nmin)")
    msg && println(", Nuctp = $(def.pos.Nuctp), Nb = $(def.pos.Nb), N = $(def.pos.N), nodes = $(def.pos.nodes), Rmax = $(grid.r[grid.N])\n")
        
    if test > 0 Z = OUTSCH!(Z, E, grid, def, adams) end
    if test > 1 Z = adams_moulton_outward!(Z, def, adams) end
    if test > 2 Z = INSCH_WKB!(Z, E, grid, def) end
    if test > 3 (ΔQ, Z) = adams_moulton_inward!(Z, def, adams) end
    if test > 4 (ΔE, Z) = adams_moulton_normalize!(Z, ΔQ, grid, def) end

    return def, Z

end
