# SPDX-License-Identifier: MIT

# author: Jook Walraven - 14-2-2023

# ==============================================================================
#                               hydrogen.jl
# ==============================================================================

# ..............................................................................
function _hydrogenic_norm(n::Int, ℓ::Int)

    o = big(1)
    n = big(n)
    ℓ = big(ℓ)

    for i=(n-ℓ):(n+ℓ)
        o *= i
    end

    o *= 2n

    return o

end
# ..............................................................................

@doc raw"""
    hydrogenic_reduced_wavefunction(atom::Atom, orbit::Orbit, grid::CamiDiff.Grid{T}) where T<:Real


Analytic expression for the hydrogenic wavefunction written in the format
``Z = \tilde{χ} + i \tilde{χ}^′``, where ``\tilde{χ}_{nℓ}(ρ)`` is
the *reduced* radial wavefunction and ``\tilde{χ}^′_{nℓ}(ρ)`` its derivative,
with ``ρ`` the radial distance to the nucleus in a.u..
The expression is evaluated for a given [`Atom`](@ref) in a given
[`Orbit`](@ref) on a given [`CamiDiff.Grid`](@extref). The argument [`Def`](@ref)
completes the definition of the problem.
```math
    \tilde{\chi}_{nl}(\rho)
    =\mathcal{N}_{nl}^{-1/2}(2Z/n)^{l+3/2}\rho^{l+1}e^{-Zρ/n}
    L_{n-l-1}^{2l+1}(2Z\rho/n)
```
where ``L_{n-l-1}^{2l+1}(2Z\rho/n)`` is the generalized Laguerre polynomial
`CamiMath.generalized_laguerreL` and
```math
    \mathcal{N}_{nl}
    = {\displaystyle \int\nolimits _{0}^{\infty}}x^{2l+2}e^{-x}
    \left[L_{n-l-1}^{2l+1}(x)\right]^{2}dx
    = \frac{2n\Gamma(n+l+1)}{\Gamma(n-l)}
```
is the norm of the wavefunction.
#### Example:
```
julia> codata = castCodata(2022);

julia> atom = castAtom(;Z=1, A=1, Q=0, msg=false);

julia> orbit = castOrbit(n=25, ℓ=10);

julia> grid = autoGrid(atom, orbit, Float64; msg=true);
CamiDiff.Grid created: exponential, Float64, rmax = 3651.58 a.u., N = 1320, h = 0.0075815, r0 = 0.164537

julia> Z = hydrogenic_reduced_wavefunction(atom, orbit, grid);
 IOP capture at generalized_laguerre_polynom(35, 21): output converted to BigInt

julia> def = castDef(grid, atom, orbit, codata; msg=true);
Def created for ¹H:25n on exponential grid of 1320 points

plot_wavefunction(Z, 1:grid.N, grid, def)
```
The plot is made using `CairomMakie`.
NB.: `plot_wavefunction` is not included in the `CamiXon` package.
![Image](../assets/H1_25n.svg)
"""
function hydrogenic_reduced_wavefunction(atom::Atom, orbit::Orbit, grid::CamiDiff.Grid{T}) where T<:Real

    N = grid.N
    n = orbit.n
    ℓ = orbit.ℓ
    r = grid.r
    Z = atom.Z

    norm = _hydrogenic_norm(n, ℓ)

    polynom = float(CamiMath.generalized_laguerre_polynom(n-ℓ-1, 2ℓ+1))

    mytype = r[N]^(ℓ+1) == Inf ? BigFloat : Float64
    a = 2Z//n
    b = mytype(a//norm)
    b = sqrt(b)
    a = mytype(a) 

    ar = a .* r
    arl = [ar[i]^ℓ for i=1:N]

    P = b .* [arl[i] * a*r[i] * exp(-Z//n*r[i]) * CamiMath.polynomial(polynom, ar[i]) for i=1:N]
    Q = b .* [arl[i] * a * exp(-Z//n*r[i]) * (((ℓ+1)-Z//n*r[i]) * CamiMath.polynomial(polynom, ar[i]) + a*r[i]*CamiMath.polynomial(polynom, ar[i]; deriv=1)) for i=1:N]

    P[1] = T(0)
    return T.(P) + im * T.(Q)

end

# ======================== bohrformula(atom, term) =============================

@doc raw"""
    bohrformula(Z::Int, n::Int)

Hydrogenic energy (in Hartree a.u.) for *atom* with *atomic number* `Z` and
*principal quantum number* `n`.
```math
    E_n = - \frac{Z^2}{2n^2}
```
#### Example:
```
julia> Z = 2; n = 4;

julia> bohrformula(Z,n)
-1//8
```
"""
bohrformula(Z::Int, n::Int) = -(1//2)*(Z//n)^2

# =================== reduce_wavefunction(Z, grid) ============================

@doc raw"""
    reduce_wavefunction(Z::Vector{Complex{T}}, grid::CamiDiff.Grid{T}) where T<:Real

Conversion from the *ordinary* radial wavefunction ``\tilde{R}_{nl}(ρ)``
to the *reduced* radial wavefuntion
```math
    \tilde{\chi}_{nl}(ρ) = ρ \tilde{R}_{nl}(ρ).
```
where ``ρ`` is the radial distance to the nucleus in a.u..
#### Example:
```
julia> atom = castAtom(Z=1, A=1, Q=0; msg=false);
julia> orbit = castOrbit(n=1, ℓ=0; msg=false);
julia> grid = autoGrid(atom, orbit, Float64);
julia> RH1s_example = [RH1s(atom.Z, grid.r[n]) for n=1:grid.N];
julia> ZH1s_example = reduce_wavefunction(RH1s_example, grid);
julia> ZH1s_generic = hydrogenic_reduced_wavefunction(atom, orbit, grid);
julia> @test ZH1s_example ≈ ZH1s_generic
Test Passed

julia> f1 = real(ZH1s_example);
julia> f2 = real(ZH1s_generic);
julia> compare_functions(f1, f2, 1:grid.N, grid)
```
The plot is made using `CairomMakie`.
NB.: `compare_functions` is not included in the `CamiXon` package.
![Image](../assets/compareXH1s.png)
```
"""
function reduce_wavefunction(Z::Vector{Complex{T}}, grid::CamiDiff.Grid{T}) where T<:Real

    R = real(Z)
    R′= imag(Z)
    r = grid.r

    χ = r .* R
    χ′= r .* R′ + R

    return χ + im * χ′

end

# =================== restore_wavefunction(Z, grid) ============================

@doc raw"""
    restore_wavefunction(Z::Vector{Complex{T}}, atom::Atom, orbit::Orbit, grid::CamiDiff.Grid{T}) where T<:Real

Conversion from the *reduced* radial wavefunction ``\tilde{\chi}_{nl}(ρ)``
to the *ordinary* radial wavefuntion ``\tilde{R}_{nl}(ρ)``,
```math
    \tilde{R}_{nl}(ρ)=\tilde{\chi}_{nl}(ρ)/ρ,
```
where ``ρ`` is the radial distance to the nucleus in a.u..
#### Example:
```
julia> atom = castAtom(Z=1, A=1, Q=0; msg=false);
julia> orbit = castOrbit(n=1, ℓ=0; msg=false);
julia> grid = autoGrid(atom, orbit, Float64);
julia> RH1s_example = [RH1s(atom.Z, grid.r[n]) for n=1:grid.N];
julia> ZH1s_example = reduce_wavefunction(RH1s_example, grid);
julia> RH1s_generic = restore_wavefunction(ZH1s_generic, atom, orbit, grid);  

julia> @test RH1s_example ≈ RH1s_generic 
Test Passed
"""
function restore_wavefunction(Z::Vector{Complex{T}}, atom::Atom, orbit::Orbit, grid::CamiDiff.Grid{T}) where T<:Real

    χ = real(Z)
    χ′= imag(Z)
    r = grid.r
    k = grid.k
    ℓ = orbit.ℓ
    Zval = atom.Z          # do not confuse the nuclear charge number Z with with the reduced wavefuntion Z

    R = χ ./ r
    R′= (χ′ .- χ ./ r) ./ r

    # obsolete: R′[1] = CamiDiff.fdiff_interpolation(R′[2:end], 0; k)   #

    if !isfinite(R′[1])  # extrapolate to r=0 to handle "division by zero"
        α = CamiDiff.fdiff_interpolation_expansion_polynom(1, k, fwd)
        Fk = CamiDiff.fdiff_expansion_weights(α, fwd, reg)
        R′[1] = LinearAlgebra.dot(Fk, R′[2:2+k])
    end
    R[1] = ℓ > 0 ? T(0) : R′[1] * convert(T, -(ℓ + 1)//Zval )  # construct R[1] from R′[1] - see OUTSCH.jl

    return R + im * R′

end

# =============================== RHqs(Z, r) ===================================

@doc raw"""
    RH1s(Z::Int, r::T) where T <:Real

Analytic expression for the *hydrogenic* 1s radial wavefunction
and its derivative in the format ``Z = \tilde{R} + i \tilde{R}^′``, where
```math
    \tilde{R}_{1s}(ρ) = Z^{3/2} 2 e^{-Zρ}
```
is the radial wavefunction and
```math
    \tilde{R}^′_{1s}(ρ) = -Z^{5/2} 2 e^{-Zρ}
```
its derivative, with ``ρ``
the radial distance to the nucleus in a.u..
#### Example:
```
atom = castAtom(Z=1, A=1, Q=0; msg=false);
orbit = castOrbit(n=1, ℓ=0; msg=false);
grid = autoGrid(atom, orbit, Float64; Nboost=1, msg=false);
def = castDef(grid, atom, orbit, codata);

RH1s_example = [RH1s(atom.Z, grid.r[n]) for n=1:grid.N];

plot_wavefunction(RH1s_example, 1:grid.N, grid, def; reduced=false)
```
The plot is made using `CairomMakie`.
NB.: `plot_function` is not included in the `CamiXon` package.
![Image](../assets/RH1s.png)
"""
function RH1s(Z::Int, r::T) where T <:Real

    one = T(1)
    two = T(2)
    Z = T(Z)

    p = two * Z^(3/2) * exp(-Z*r)
    q = -two * Z^(5/2) * exp(-Z*r)

    return p + im * q

end

# =============================== RH2p(Z, r) ===================================

@doc raw"""
    RH2s(Z::Int, r::T) where T <:Real

Analytic expression for the *hydrogenic* 1s reduced radial wavefunction
and its derivative in the format ``Z = \tilde{R} + i \tilde{R}^′``, where
```math
    \tilde{R}_{2s}(ρ)=\left(Z/2\right)^{3/2}(1-Zρ/2)2e^{-Zρ/2}
```
is the radial wavefunction and
```math
    \tilde{R}_{2s}(ρ)=-\left(Z/2\right)^{5/2}(2-Zρ/2)2e^{-Zρ/2}
``` its derivative, with ``ρ`` the radial distance to the nucleus in a.u..
#### Example:
```
atom = castAtom(Z=1, A=1, Q=0; msg=false);
orbit = castOrbit(n=2, ℓ=0; msg=false);
grid = autoGrid(atom, orbit, Float64; Nboost=1, msg=false);
def = castDef(grid, atom, orbit, codata; msg=false);

RH2s_example = [RH2s(atom.Z, grid.r[n]) for n=1:grid.N];

plot_wavefunction(RH2s_example, 1:grid.N, grid, def; reduced=false)
```
The plot is made using `CairomMakie`.
NB.: `plot_wavefunction` is not included in the `CamiXon` package.
![Image](../assets/RH2s.png)
"""
function RH2s(Z::Int, r::T) where T <:Real

    one = T(1)
    two = T(2)
    Z = T(Z)

    P = two * (Z/two)^(3/2) * (one-Z*r/two) * exp(-Z*r/two)
    Q = -two * (Z/two)^(5/2) * (two - Z*r/two) * exp(-Z*r/two)

    return P + im * Q

end

# =============================== RH2p(Z, r) ===================================

@doc raw"""
    RH2p(Z::Int, r::T) where T <:Real

Analytic expression for the *hydrogenic* 1s reduced radial wavefunction
and its derivative in the format ``Z = \tilde{R} + i \tilde{R}^′``, where
```math
    \tilde{R}_{2p}(ρ)=\left(Z/2\right)^{3/2}\sqrt{1/3}(Zρ/2)2e^{-Zρ/2}
```
is the radial wavefunction and
```math
    \tilde{R}_{2p}(ρ)=\left(Z/2\right)^{3/2}\sqrt{1/3}(1-Zρ/2)2e^{-Zρ/2}
``` its derivative, with ``ρ`` the radial distance to the nucleus in a.u..
#### Example:
```
atom = castAtom(Z=1, A=1, Q=0; msg=false);
orbit = castOrbit(n=2, ℓ=1; msg=false);
grid = autoGrid(atom, orbit, Float64; Nboost=1, msg=false);
def = castDef(grid, atom, orbit, codata);

RH2p_example = [RH2p(atom.Z, grid.r[n]) for n=1:grid.N];

plot_wavefunction(RH2p_example, 1:grid.N, grid, def; reduced=false)
```
The plot is made using `CairomMakie`.
NB.: `plot_wavefunction` is not included in the `CamiXon` package.
![Image](../assets/RH2p.png)
"""
function RH2p(Z::Int, r::T) where T <:Real

    one = T(1)
    two = T(2)
    tri = T(3)
    Z = T(Z)
   
    P = two * (Z/two)^(3/2) * sqrt(1/tri) * (Z*r/two) * exp(-Z*r/two)
    Q = two * (Z/two)^(5/2) * sqrt(1/tri) * (one - Z*r/two) * exp(-Z*r/two)

    return P + im * Q

end

#χHe1s(r) = 4.0 * sqrt(2) * r * exp(-2.0r) + im * 4.0 * sqrt(2) * exp(-2.0r) * (1 - 2.0r)
#gridHe1s(grid) = [χHe1s(grid.r[n]) for n=1:grid.N]

@doc raw"""
    silvera_goldman_triplet(r::T) where T<:Real

Parametrization in *Hartree* a.u. of the triplet ``(^{3}\Sigma_{u}^{+})`` potential of the electronic 
ground state of ``\mathrm{H}_{2}`` (Eh = 219474.6 cm-1), 

```math
   V_{t}(r)=\mathrm{exp}\left(0.09678-1.10173\thinspace r-0.0394\thinspace r^{2}\right)+F(r)\left(-6.5\thinspace r^{-6}-124\thinspace r^{-8}-3285r^{-10}\thinspace \right)
```
```math
   \mathrm{where}\ \ \ \ \ \ \ \  F(r) = \begin{cases}
\mathrm{exp}\left[-\left(\frac{10.04}{r}-1\right)^{2}\right] & \mathrm{for}\,\,\,r<10.04\,\mathrm{a.u.}\\
1 & \mathrm{for}\,\,\,r<10.04\,\mathrm{a.u.}
\end{cases}
```
see I.F. Silvera, - Rev. Mod. Phys., 52, 393 (1980).
"""
function silvera_goldman_triplet(r::T) where T<:Real 
# ====================================================== 
#    I.F. Silvera, - Rev. Mod. Phys., 52, 393 (1980)
#
#    Hartree: 219474.6 cm-1    
# ======================================================

    Fr = r > 10.04 ? 1.0 : exp(-(10.04/r-1.0)^2)
   
    o = exp(0.09678 - 1.10173r - 0.03945r^2) + Fr * (-6.5/r^6-124.0/r^8-3285.0/r^10)

    return o

end

@doc raw"""
    silvera_goldman_exchange(r::T) where T<:Real

Parametrization in *Hartree* a.u. of the exchange energy difference between singlet ``(^{3}\Sigma_{u}^{+})`` 
and triplet ``(^{1}\Sigma_{g}^{+})`` potentials  of ``\mathrm{H}_{2}`` (Eh = 219474.6 cm-1), 
```math
    J(r)=\mathrm{exp}\left(-0.288-0.275\thinspace r-0.176\thinspace r^{2}+0.0068\thinspace r^{3}\right)
```
see I.F. Silvera, - Rev. Mod. Phys., 52, 393 (1980).
"""
function silvera_goldman_exchange(r::T) where T<:Real 
# ====================================================== 
#    I.F. Silvera, - Rev. Mod. Phys., 52, 393 (1980)
#
#    Hartree: 219474.6 cm-1    
# ====================================================== 
   
    o = r < 18.0 ? exp(-0.288 - 0.275 * r - 0.176 * r^2 + 0.0068 * r^3) - 1.5242e-10 : 0.0

    return o

end

@doc raw"""
    silvera_goldman_singlet(r::T) where T<:Real 

Parametrization in *Hartree* a.u. of the singlet ``(^{1}\Sigma_{g}^{+})`` potential of ``\mathrm{H}_{2}`` (Eh = 219474.6 cm-1), 
```math
   V_{s}(r)=V_{t}(r)-J(r)
```
where ``V_{t}`` is the triplet potential ([`silvera_goldman_triplet`](@ref)) and ``J(r)`` is the exchange contribution ([`silvera_goldman_triplet`](@ref)).
"""
function silvera_goldman_singlet(r::T) where T<:Real 

    U = silvera_goldman_triplet(r)
    J = silvera_goldman_exchange(r)
    o = U - J

    return o

end


@doc raw"""
    silvera_goldman_potential(grid::CamiDiff.Grid{T}; S=0) where T<:Real

CamiDiff.Grid representation in *Hartree* a.u. of the singlet (S=0) and triplet (S=1) potentials of ``\mathrm{H}_{2}``,
```math
    \mathcal{V}(r)=V_{D}(r)+J(r)\mathbf{s}_{1}\cdot\mathbf{s}_{2},
```
where ``\mathbf{S} = \mathbf{s}_{1}+\mathbf{s}_{2}`` and 
```math
    V_{D}(r)=\frac{1}{4}[V_{s}(r)+3V_{t}(r)] \ \mathrm{and}\ J(r)=V_{t}(r)-V_{s}(r)
```
are known as the direct and exchange contribution, respectively; ``V_{s}:`` see [`silvera_goldman_singlet`](@ref), 
``V_{t}:`` see [`silvera_goldman_triplet`](@ref), ``J(r):`` see [`silvera_goldman_exchange`](@ref).
    
see I.F. Silvera, - Rev. Mod. Phys., 52, 393 (1980).
"""
function silvera_goldman_potential(grid::CamiDiff.Grid{T}; S=0) where T<:Real

    if iszero(S)
        o = [silvera_goldman_singlet(grid.r[i]) for i=1:grid.N]    # singlet potential
    else
        o = [silvera_goldman_triplet(grid.r[i]) for i=1:grid.N]    # triplet potential
    end

    return o
    
end

@doc raw"""
    rotbarrier(grid::CamiDiff.Grid{T}; ℓ=0) where T<:Real

CamiDiff.Grid representation of rotational barrier potential in wavenumber notation,
```math
V_{rot}(r) = \frac{\ell(\ell+1)}{r^{2}},
``` 
where ℓ is the rotational quantum number.
"""
function rotbarrier(grid::CamiDiff.Grid{T}; ℓ=0) where T<:Real

    if ℓ > 0
        num = convert(T, ℓ*(ℓ + 1))
        o = [num*(grid.r[i])^-2 for i=1:grid.N]
    else
        o = zeros(T, grid.N)
    end

    return o
    
end