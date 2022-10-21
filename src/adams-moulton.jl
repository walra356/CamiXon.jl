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

    norm = abs(real(Z[Nuctp]))

    Z[1:Nuctp] /= norm    # set amplitude at u.c.t.p. to +1/-1 (nodes even/odd)

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
    Nb = def.pos.Nb

    for n=Nb-1:-1:Nuctp
        _prepend!(Z2, n, adams.Minv, adams.G, def.am, k)
    end

    norm = abs(real(Z2[1]))
    Z2 /= norm
    Z2 *= sgn             # set amplitude at u.c.t.p. to +1/-1 (nodes even/odd)

    ΔQ = imag(Z[Nuctp]) - imag(Z2[1])

    Z[Nuctp:N] = Z2

    def.pos.Na = get_Na(Z, def)
    def.pos.Nb = get_Nb(Z, def)

    return ΔQ, Z

end
function adams_moulton_inward_old(E::T, grid::Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real

    N = grid.N
    k = grid.k
    Z = adams.Z
    Nuctp = def.pos.Nuctp
    Nb = def.pos.Nb
    sgn = isodd(def.orbit.n′) ? -1.0 : 1.0

    Z2 = INSCH(E, grid, def, adams)

    for n=Nb-1:-1:Nuctp
        _prepend!(Z2, n, adams.Minv, adams.G, def.am, k)
    end

    norm = abs(real(Z2[1]))
    Z2 /= norm
    Z2 *= sgn             # set amplitude at u.c.t.p. to +1/-1 (nodes even/odd)

    ΔQ = imag(Z[Nuctp]) - imag(Z2[1])

    Z[Nuctp:N] = Z2

    def.pos.Na = get_Na(Z, def)
    def.pos.Nb = get_Nb(Z, def)

    return ΔQ, Z

end
function adams_moulton_inward_WJ(E::T, grid::Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real

    N = grid.N
    k = grid.k
    Z = adams.Z
    Nuctp = def.pos.Nuctp
    sgn = isodd(def.orbit.n′) ? -1.0 : 1.0

    Z2 = INSCH(E, grid, def, adams)

    for n=N-k-1:-1:Nuctp
        _prepend!(Z2, n, adams.Minv, adams.G, def.am, k)
    end

    norm = abs(real(Z2[1]))
    Z2 /= norm
    Z2 *= sgn             # set amplitude at u.c.t.p. to +1/-1 (nodes even/odd)

    ΔQ = imag(Z[Nuctp]) - imag(Z2[1])

    Z[Nuctp:N] = Z2

    def.pos.Na = get_Na(Z, def)
    def.pos.Nb = get_Nb(Z, def)

    return ΔQ, Z

end

# ======================= adams_moulton_normalized =============================

@doc raw"""
    adams_moulton_normalized(Z::Vector{Complex{T}}, ΔQ::T, grid::Grid{T}, def::Def{T}) where T<:Real

"""
function adams_moulton_normalized(Z::Vector{Complex{T}}, ΔQ::T, grid::Grid{T}, def::Def{T}) where T<:Real

    Nuctp = def.pos.Nuctp
    #Na = def.pos.Na
    #Nb = def.pos.Nb
    #k = def.k

    #norm = grid_integration(real(Z) .^2, Na-k, Nb+k, grid)

    norm = grid_integration(real(Z) .^2, 1, N, grid)

    ΔE = ΔQ * abs(real(Z[Nuctp])) / T(2)

    return ΔE/norm, Z/sqrt(norm)

end

# =============== adams_moulton_patch(E, grid, def, adams) =====================

@doc raw"""
    adams_moulton_patch(Z::Vector{Complex{T}}, def::Def{T}, adams::Adams{T}) where T<:Real

"""
function adams_moulton_patch(Z::Vector{Complex{T}}, def::Def{T}, adams::Adams{T}) where T<:Real

    k = def.k

    Z2 = copy(Z[2k+1:3k])

    for n=2k:-1:1
        _prepend!(Z2, n, adams.Minv, adams.G, def.am, k)
    end

    Z[1:3k] = Z2

    P = fdiff_interpolation(real(Z[2:end]),0; k=4)
    Q = fdiff_interpolation(imag(Z[2:end]),0; k=4)

    Z[1] = P + im * Q

    return Z

end

# =============== adams_moulton_solve(E, grid, def, adams) =====================

@doc raw"""
    adams_moulton_solve(E::T, grid::Grid{T}, def::Def{T}, adams::Adams) where T<:Real

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
function adams_moulton_solve(E::T, grid::Grid{T}, def::Def{T}, adams::Adams) where T<:Real

    adams = updateAdams!(adams, E, grid, def)
        Z = adams_moulton_outward(def, adams)
    ΔQ, Z = adams_moulton_inward(E, grid, def, adams)
    ΔE, Z = adams_moulton_normalized(Z, ΔQ, grid, def)

    if def.pos.Na == def.k + 1
         Z = adams_moulton_patch(Z, def, adams)
    end

    for n ∈ eachindex(Z)
        adams.Z[n] = Z[n]
    end

    return adams, ΔE, Z

end

# ======================= adams_moulton_master sector ==========================

# ..............................................................................
function _set_bounds!(init::NTuple{4,T}, grid::Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real

    n′= def.orbit.n′     # radial quantum number (number of nodes)

    (Emin, E, Emax, ΔE) = init

    Z = adams.Z

    nodes = def.pos.nodes

    i = 0
    while nodes ≠ n′
        if nodes < n′
            Emin = max(E,Emin)
            E = T(0.9E)
        else
            Emax = min(E,Emax)
            E = T(1.1E)
        end
        adams, ΔE, Z = adams_moulton_solve(E, grid, def, adams)
        nodes = def.pos.nodes = count_nodes(Z, def)
        i += 1
    end

    init = (Emin, E, Emax, ΔE)

    return i, def, adams, init, Z

end
# ..............................................................................
function _message(i::Int, imax::Int, init::NTuple{4,T}, def::Def{T}; modus="prepare") where T<:Real

    nodes = def.pos.nodes
       n′ = def.orbit.n′
    codata = def.codata

    nodes == n′ || error("Error: nodes condition violated (nodes = $(nodes) ≠  $(n′)")

    (Emin, E, Emax, ΔE) = init

    Nuctp = def.pos.Nuctp
    Ntot = def.pos.N
    cnts = "iterations"
    f = convertUnit(abs(ΔE), codata) # default input (Hartree) and output (xHz)
    strHz = strValue(f)
    imax = modus == "prepare" ? "" : i > imax-1 ? "(maximum value)" : ""
    strΔE = @sprintf "ΔE = %.17g %s" abs(ΔE) " Hartree (" * strHz *") - absolute convergence error\n"
    strΔErel = @sprintf "ΔE/E = %.17g %s" abs(ΔE/E) " - relative convergence error\n"

    msg = "\n" * modus * "_adams_moulton (" * string(def.T) * "): \nnode condition satified after "
    msg *= "$i " * cnts * imax
    msg *= ": nodes = $(nodes), "
    msg *= "Nuctp = $(Nuctp), "
    msg *= "Ntot = $(Ntot) \n"
    msg *= @sprintf "E = %.17g %s\n" E "Hartree"
    msg *= strΔE
    msg *= strΔErel

end
# ..............................................................................
function _strΔt(tstop::T, tstart::T) where T<:Real

    Δt = tstop-tstart
    str = Δt > 1.0  ? (repr(Δt, context=:compact => true) * " sec")      :
          Δt > 1e-3 ? (repr(Δt*1e3, context=:compact => true) * " msec") :
                      (repr(Δt*1e6, context=:compact => true) * " μsec")
    return str

end
# ..............................................................................

# ============================ adams_moulton_prepare ===========================

@doc raw"""
    adams_moulton_prepare(E::T, grid::Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real

Solves the Schrödinger equation for an atom defined by `def` for energy `E`
on grid the `grid` with the Adams-Moulton method defined by `adams`. `E` is
adjusted until the wavefunction has the correct number of `n′` nodes.

#### Example:

```
Ecal, grid, def, adams = demo_hydrogen(n=1, ℓ=0);
    Def created for hydrogen 1s on exponential grid of 100 points

E = 1.5Ecal
msg, adams, init, Z = adams_moulton_prepare(E, grid, def, adams);
    Ecal = -0.5; E = -0.75; 0 nodes

plot_wavefunction(Z, 1:def.pos.N, grid, def; reduced=false)
```
The plot is made using `CairomMakie`. Note the discontinuity in the derivative.
NB.: `plot_wavefunction` is not included in the `CamiXon` package.

![Image](./assets/hydrogen-1s-prepared.png)
"""
function adams_moulton_prepare(E::T, grid::Grid{T}, def::Def{T}, adams::Adams{T}) where T<:Real

    n′= def.orbit.n′
    N = def.pos.N
    v = def.pot
    s = def.scr

    pot = v .+ s

    Emax = pot[N]
    Emin = minimum(pot[2:N])

    adams, ΔE, Z = adams_moulton_solve(E, grid, def, adams)

    nodes = def.pos.nodes = count_nodes(Z, def)
    nodes > n′ ? Emax = E : nodes < n′ ? Emin = E : false

    init = (Emin, E, Emax, ΔE)

    i, def, adams, init, Z = _set_bounds!(init, grid, def, adams)

    msg = _message(i, 0, init, def; modus="prepare")

    return msg, adams, init, Z

end

# ========================= adams_moulton_iterate ==============================

@doc raw"""
    adams_moulton_iterate(init::NTuple{4,T}, grid::Grid{T}, def::Def{T}, adams::Adams{T}; imax=25, Δν=Value(1,"kHz")) where T<:Real

Solves the Schrödinger equation for an atom defined by `def` for energy `E`
on grid the `grid` with the Adams-Moulton method defined by `adams`; `E` is
adjusted in an iteration procedure until convergence is reached within the
convergence goal `Δν` is reached (limited to a maximum of `imax` iterations).

#### Example:
```
Ecal, grid, def, adams = demo_hydrogen(n=1, ℓ=0);
    Def created for hydrogen 1s on exponential grid of 100 points

E = 1.5Ecal;
msg1, adams, init, Z = adams_moulton_prepare(E, grid, def, adams);
println("Ecal = $Ecal; E = $(init[2]); $(def.pos.nodes) nodes")
    Ecal = -0.5; E = -0.75; 0 nodes

msg2, adams, init, Z = adams_moulton_iterate(init, grid, def, adams; Δν=Value(1,"MHz"), imax=25)
println("Ecal = $Ecal; E = $(init[2]); $(def.pos.nodes) nodes")
    Ecal = -0.5; E = -0.49999997841850014; 0 nodes

plot_wavefunction(Z, 1:def.pos.N, grid, def; reduced=false)
```
The plot is made using `CairomMakie`.
NB.: `plot_wavefunction` is not included in the `CamiXon` package.
![Image](./assets/hydrogen-1s.png)
"""
function adams_moulton_iterate(init::NTuple{4,T}, grid::Grid{T}, def::Def{T}, adams::Adams{T}; imax=25, Δν=Value(1,"kHz")) where T<:Real

    n′= def.orbit.n′  # radial quantum number (number of nodes)
    codata = def.codata

    (Emin, E, Emax, ΔE) = init

    test = convertUnit(Δν.val, codata; unitIn=Δν.unit, unitOut="Hartree")

    test = T == BigFloat ? convert(T,test.val) : convert(T,test.val)

    Z = adams.Z

    i=0
    while abs(ΔE) > test
        adams, ΔE′, Z = adams_moulton_solve(E, grid, def, adams)
        nodes = def.pos.nodes = count_nodes(Z, def)
        if nodes > n′
            Emax = min(E,Emax)
            E = (E-ΔE+Emax)/2
        elseif nodes < n′
            Emin = max(E,Emin)
            E = (E-ΔE+Emin)/2
        else
            ΔE=ΔE′
            E = iseven(nodes) ? E+ΔE : E-ΔE
        end
        i = i < imax ? i+1 : break
    end

    init = (Emin, E, Emax, ΔE)

    def.pos.nodes = count_nodes(Z, def)

    msg = _message(i, imax, init, def; modus="iterate")

    return msg, adams, init, Z

end

# ===================== adams_moulton_master ===================================

@doc raw"""
    adams_moulton_master(E, grid, def, adams; Δν=Value(1,"kHz"), imax=25, msg=true)

Solves the Schrödinger equation for an atom defined by `def` for energy `E`
on grid the `grid` with the Adams-Moulton method defined by `adams`.

`Δν`: convergence goal

`imax`: maximum number of iterations

#### Example:
```
Ecal, grid, def, adams = demo_hydrogen(n=1, ℓ=0);
    Def created for hydrogen 1s on exponential grid of 100 points

E = 1.5Ecal;
E, def, adams, Z = adams_moulton_master(E, grid, def, adams; Δν=Value(1,"kHz"), imax=25, msg=true);
plot_wavefunction(Z, 1:def.pos.N, grid, def; reduced=false)
```
The plot is made using `CairomMakie`.
NB.: `plot_wavefunction` is not included in the `CamiXon` package.
![Image](./assets/hydrogen-1s.png)
"""
function adams_moulton_master(E, grid, def, adams; Δν=Value(1,"kHz"), imax=25, msg=true)

    t1 = time()
    msg && println("\nt = " * _strΔt(t1,t1))
    msg1, adams, init, Z = adams_moulton_prepare(E, grid, def, adams)
    t2 = time()
    msg && println(msg1 * "\n" * "preparation time: " * _strΔt(t2,t1))
    msg && println("\nt = " * _strΔt(t2,t1))
    msg2, adams, init, Z = adams_moulton_iterate(init, grid, def, adams; Δν, imax)
    t3 = time()
    msg && println(msg2 * "\n" * "iteration time: " * _strΔt(t3,t2))
    msg && println("\nt = " * _strΔt(t3,t1))

    (Emin, E, Emax, ΔE) = init

    msg && println("\ngrid special points: Na = $(def.pos.Na), Nlctp = $(def.pos.Nlctp), Nuctp = $(def.pos.Nuctp), Nb = $(def.pos.Nb), N = $(def.pos.N)")

    def.orbit.n′ == count_nodes(Z, def) || error("Error: node condition violated - n′ = $(def.orbit.n′), nodes = $(nodes)")

    return E, def, adams, Z

end
