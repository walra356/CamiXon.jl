# SPDX-License-Identifier: MIT

# Copyright (c) 2025 Jook Walraven <69215586+walra356@users.noreply.github.com> and contributors

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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
    ℓ = def.spinorbit.ℓ
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
    
    Nb+k ≤ N || error("Error: grid boost required (increase N)")

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

    Q1 = imag(Z[Nuctp])
    
    ΔQ = Q1 - Q0

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
    two = T(2)

    norm = CamiDiff.grid_integration(real(Z) .^2, grid, 1, N)

    ΔE = ΔQ * abs(real(Z[Nuctp]))

    ΔE = ΔE/norm/two

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
![Image](../../assets/hydrogen-1s-prepared.png)
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

    n = def.spinorbit.n      # principal quantum number
    n′= def.spinorbit.n′     # radial quantum number (number of nodes)
    Zc= def.atom.Zc      # Rydberg charge

    E = iszero(scr[1]) ? T(-(Zc//n)^2//2) : T(E)
    
    for n ∈ eachindex(def.pot)
        def.scr[n] = scr[n]
        def.potscr[n] = def.pot[n] + scr[n]
    end
    
    inE = castInE(E, def)     # inE =(E[Nlctp], E[Nmin], E[Nuctp], ΔE=0.0)

    adams = castAdams(inE.E, grid, def)
    
    Z = zeros(Complex{grid.T}, grid.N)  

    adams, ΔE, Z = adams_moulton_solve!(Z, inE.E, grid, def, adams)
    nodes = def.pos.nodes
    if msg
        str =  Printf.@sprintf "%.20g, " inE.Emin 
        str *= Printf.@sprintf "%.20g, " inE.E
        str *= Printf.@sprintf "%.20g, " inE.Emax
        str *= Printf.@sprintf "%.4g" inE.ΔE
        println("\n--- adams_moulton_nodes!:\nstart solution: $(nodes) nodes - inE = (" * str * ")")
        print("start node search:    count = ")
    end
    
    i = 0
    while def.pos.nodes ≠ n′
        msg && print(def.pos.nodes, ",")
        adams, ΔE, Z = adams_moulton_solve!(Z, inE.E, grid, def, adams) # inE.ΔE not used for nodes ≠ n′
        inE = inE!(inE, T(0), def)
        i += 1
        i < imax || break
    end

    Q = imag(Z)
    E = inE.E

    while sign(Q[def.pos.Nuctp-2]) ≠ sign(Q[def.pos.Nuctp+2])
        E *= T(0.9)
        adams, ΔE, Z = adams_moulton_solve!(Z, E, grid, def, adams)
        Q = imag(Z)
        inE.E = def.pos.nodes == n′ ? E : inE.E/T(0.9)
        i += 1
        i < imax || break
    end

    adams, ΔE, Z = adams_moulton_solve!(Z, inE.E, grid, def, adams)

    nodes = def.pos.nodes
    
    nodes == n′ || println("nodes failed")   

    inE.ΔE = ΔE  # ΔE from here on in use for reporting only (nodes == n′)
    
    if msg
        str1 = nodes == n′ ? "$(nodes) - node search completed\n" : error("Error: after $i iterations nodes = $(nodes) - should be $(n′) (increase imax?)")
        str =  Printf.@sprintf "%.20g, " inE.Emin 
        str *= Printf.@sprintf "%.20g, " inE.E
        str *= Printf.@sprintf "%.20g, " inE.Emax
        str *= Printf.@sprintf "%.4g" inE.ΔE
        println(str1 * "initiate ΔE: $(nodes) nodes - inE = (" * str * ")")
    end
    
    t2 = time()
    adams_moulton_report_nodes(i, inE, grid, def, _strΔt(t2,t1); unitIn="Hartree", msg)

    return def, adams, inE, Z

end

# --------------------------------------------------------------------------------------------------------------
#            adams_moulton_iterate!(Z, inE, grid, def, adams; imax=25, ϵ=1e-6, msg=true)
# --------------------------------------------------------------------------------------------------------------

@doc raw"""
    adams_moulton_iterate!(Z::Vector{Complex{T}}, inE::InE{T}, grid::CamiDiff.Grid{T}, def::Def{T}, adams::Adams{T}; imax=25, ϵ=1e-6, msg=true) where T<:Real
    
"""
function adams_moulton_iterate!(Z::Vector{Complex{T}}, inE::InE{T}, grid::CamiDiff.Grid{T}, def::Def{T}, adams::Adams{T}; imax=25, ϵ=1.0e-6, msg=true) where T<:Real
    
    t1 = time()

    n′= def.spinorbit.n′     # radial quantum number (number of nodes)
    nodes = def.pos.nodes
    #two = T(2)
    #pot = def.potscr

    if msg
        str =  Printf.@sprintf "%.20g, " inE.Emin 
        str *= Printf.@sprintf "%.20g, " inE.E
        str *= Printf.@sprintf "%.20g, " inE.Emax
        str *= Printf.@sprintf "%.4g" inE.ΔE
        println("\n--- adams_moulton_iterate!:\nresume solution: $(nodes) nodes - inE = (" * str * ")")
    end
    
    Emin = inE.Emin
    E = inE.E
    Emax = inE.Emax
    ΔE = 1.1e-3 * E

   (Emin ≤ E ≤ Emax) || error("Error: Emin ≤ E ≤ Emax interval violated")     
 
    i = 0
    while abs(ΔE/E) > 1e-3 # convergence goal
        adams, ΔE, Z = adams_moulton_solve!(Z, E, grid, def, adams)
        inE = inE!(inE, ΔE, def)
        E = inE.E
        i += 1
        i < imax || break
    end
    abs(ΔE/E) ≤ 1e-3 || error("Error: after $i iterations ΔE/E = $(ΔE/E) - should be < 1e-3?)")
    while abs(big(ΔE)/big(E)) > ϵ # convergence goal
        ref = abs(ΔE)
        adams, ΔE, Z = adams_moulton_solve_refine!(Z, E, grid, def, adams)
        inE = inE!(inE, ΔE, def)
        E = inE.E
        i += 1
        i < imax || break
        abs(ΔE) < ref || break
    end
    
    if msg
        str =  Printf.@sprintf "%.20g, " inE.Emin 
        str *= Printf.@sprintf "%.20g, " inE.E
        str *= Printf.@sprintf "%.20g, " inE.Emax
        str *= Printf.@sprintf "%.4g" inE.ΔE
        println(                        "     upgrade ΔE: $(nodes) nodes - inE = (" * str * ")")
    end
    
    t2 = time()
    adams_moulton_report_iterate(i, imax, inE, ϵ, grid, def, _strΔt(t2,t1); unitIn="Hartree", msg)
    msg && println("grid range and special points: rmax = $(round(Float64(grid.r[grid.N]))) a.u.:  " * listPos(def.pos) )
    
    return def, adams, inE, Z

end

# --------------------------------------------------------------------------------------------------------------
#          adams_moulton_report_nodes(it::Int, inE::InE{T}, grid::CamiDiff.Grid{T}, def::Def{T}, strΔT::String; unitIn="Hartree", msg=true) where T<:Real
# --------------------------------------------------------------------------------------------------------------

@doc raw"""
    adams_moulton_report_nodes(i::Int, inE::InE{T}, grid::CamiDiff.Grid{T}, def::Def{T}, strΔT::String; unitIn="Hartree", msg=true) where T<:Real

"""
function adams_moulton_report_nodes(i::Int, inE::InE{T}, grid::CamiDiff.Grid{T}, def::Def{T}, strΔT::String; unitIn="Hartree", msg=true) where T<:Real

    ΔE = inE.ΔE
    Δf = convertUnit(abs(ΔE), def.codata)
    strΔf = " (" * strValue(Δf) * ")"
    strΔE = repr(ΔE, context=:compact => true)
    strΔErel = repr(ΔE/inE.E, context=:compact => true)
    n′= def.spinorbit.n′
    n = def.pos.nodes

    str = "\nadams_moulton_nodes: report for " * _defspecs(grid, def) * " (using $T)\n"
    str *= n == n′ ? "target of $n nodes achieved after $(i) iterations in " * strΔT * "\n" :
                     "found $n nodes in $(i) iterations - Error: $(n′) nodes expected - increase imax and/or N\n" 
    str *= Printf.@sprintf "    binding energy: E = %.20g %s \n" inE.E unitIn
    str *= ΔE ≠ 0 ? "absolute precision: ΔE = " * strΔE * " " * unitIn * strΔf * "\n" :
                    "absolute precision: ΔE = 0 (exact under $T precision)\n"
    str *= ΔE ≠ 0 ? "relative precision: ΔE/E = " * strΔErel * ""                   :
                    "relative precision: ΔE/E = 0 (exact under $T precision)"

    return msg ? println(str) : str

end

# --------------------------------------------------------------------------------------------------------------
#          adams_moulton_report_iterate(i, inE, ϵ, grid, def::Def{T}, strΔT::String; unitIn="Hartree", ϵ, msg=true)
# --------------------------------------------------------------------------------------------------------------

@doc raw"""
    adams_moulton_report_iterate(i::Int, imax::Int, inE::InE{T}, ϵ, grid::CamiDiff.Grid{T}, def::Def{T}, strΔT::String; unitIn="Hartree", msg=true) where T<:Real

"""
function adams_moulton_report_iterate(i::Int, imax::Int, inE::InE{T}, ϵ, grid::CamiDiff.Grid{T}, def::Def{T}, strΔT::String; unitIn="Hartree", msg=true) where T<:Real

    ϵ = T(ϵ)
    ΔE = inE.ΔE
    ϵv = abs(ΔE/inE.E)
    Δf = convertUnit(abs(ΔE), def.codata)
    strΔf = " (" * strValue(Δf) * ")"
    strΔE = repr(ΔE, context=:compact => true)
    strΔErel = repr(ΔE/inE.E, context=:compact => true)
    n′= def.spinorbit.n′
    n = def.pos.nodes
    strNodes = n′ == n ? "Passed node test ($n nodes); " : "Failed nodes test ($(n′) ≠ $n); "

    str = "\nadams_moulton_iterate!: report for " * _defspecs(grid, def) * " (using $T)\n"
    if ϵv < ϵ 
        str *= "reached covergence goal ($(Float64(ϵv)) < ϵ = $(Float64(ϵ))) after $(i) iterations in " * strΔT * "\n" 
    else
        str *= i == imax ? ("Warning: stopped at i=imax=$(i) after " * strΔT * " - failed to reach covergence goal ($(Float64(ϵv)) > ϵ = $(Float64(ϵ))) - increase imax and/or N (or increase ϵ)\n") : 
                           ("Warning: reached numerical resolution limit ($(Float64(ϵv))) after $(i) iterations in " * strΔT * " - consider BigFloat resolution\n")
    end
    str *= Printf.@sprintf "    binding energy: E = %.20g %s \n" inE.E unitIn
    str *= ΔE ≠ 0 ? "absolute precision: ΔE = " * strΔE * " " * unitIn * strΔf * "\n" :
                    "absolute precision: ΔE = 0 (exact under $T precision)\n"
    str *= ΔE ≠ 0 ? "relative precision: ΔE/E = " * strΔErel * ""                   :
                    "relative precision: ΔE/E = 0 (exact under $T precision)"

    return msg ? println(str) : str

end

function test_adams_moulton(E::Real, scr::Vector{T}, grid::CamiDiff.Grid{T}, def::Def{T}; test=5, msg=false) where T<:Real
    
    t1 = time()

    E = T(E)

    n′= def.spinorbit.n′     # radial quantum number (number of nodes)
    
    for n ∈ eachindex(def.pot)
        def.scr[n] = scr[n]
        def.potscr[n] = def.pot[n] + scr[n]
    end
    
    inE = castInE(E, def) # inE = (E[Nlctp], E[Nmin], E[Nuctp], ΔE=0.0)
    
    E = inE.E

    adams = castAdams(E, grid, def)
    
    Z = zeros(Complex{grid.T}, grid.N) 

    updatePos!(def.pos, E, def.potscr, grid)
    adams = updateAdams!(adams, E, grid, def)
    msg && print("grid special points: Na = $(def.pos.Na), Nlctp = $(def.pos.Nlctp), Nmin = $(def.pos.Nmin)")
    msg && println(", Nuctp = $(def.pos.Nuctp), Nb = $(def.pos.Nb), N = $(def.pos.N), nodes = $(def.pos.nodes), rmax = $(grid.r[grid.N])\n")
        
    if test > 0 Z = OUTSCH!(Z, E, grid, def, adams) end
    if test > 1 Z = adams_moulton_outward!(Z, def, adams) end
    if test > 2 Z = INSCH_WKB!(Z, E, grid, def) end
    if test > 3 (ΔQ, Z) = adams_moulton_inward!(Z, def, adams) end
    if test > 4 (ΔE, Z) = adams_moulton_normalize!(Z, ΔQ, grid, def) end

    return def, Z

end
