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
#                               CoulombIntegrals.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#           a_direct(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)
# ------------------------------------------------------------------------------

@doc raw"""
    a_direct(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)

Coulomb angular integral - direct part:

```math
a^{k}(lm_{l};l^{\prime}m_{l^{\prime}})=(-)^{m_{l}+m_{l^{\prime}}}
(2l+1)(2l^{\prime}+1)\left(\begin{array}{ccc}
l & k & l\\
0 & 0 & 0
\end{array}\right)\left(\begin{array}{ccc}
l & k & l\\
-m_{l} & 0 & m_{l}
\end{array}\right)\left(\begin{array}{ccc}
l^{\prime} & k & l^{\prime}\\
0 & 0 & 0
\end{array}\right)\left(\begin{array}{ccc}
l^{\prime} & k & l^{\prime}\\
-m_{l^{\prime}} & 0 & m_{l^{\prime}}
\end{array}\right)
```
#### Example:
```
a_direct(2,1,1,2,2)
    2//35

a_direct(6,3,2,3,-1)
    -250//20449
```
"""
function a_direct(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)

    Base.iseven(k) || return 0
    0 ≤ k ≤ 2min(l, l′) || return 0
    k > 0 || return 1

    Δ1 = CamiMath.triangle_coefficient(l, k, l)
    Δ2 = CamiMath.triangle_coefficient(l′, k, l′)
    T = CamiMath._Racah_sqrt2(l, 0, k, 0, l, 0) * CamiMath._Racah_sqrt2(l, -ml, k, 0, l, ml) * CamiMath._Racah_sqrt2(l′, 0, k, 0, l′, 0) * CamiMath._Racah_sqrt2(l′, -ml′, k, 0, l′, ml′)
    S = CamiMath._Racah_sum(l, 0, k, 0, l) * CamiMath._Racah_sum(l, -ml, k, 0, l) * CamiMath._Racah_sum(l′, 0, k, 0, l′) * CamiMath._Racah_sum(l′, -ml′, k, 0, l′)

    a = (2l + 1) * (2l′ + 1) * Δ1 * Δ2 * S
    o = a * a * T
    num = Int(sqrt(numerator(o)))
    den = Int(sqrt(denominator(o)))

    o = sign(S) * num // den

    return o

end

# ------------------------------------------------------------------------------
#           b_exchange(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)
# ------------------------------------------------------------------------------

@doc raw"""
    b_exchange(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)

Coulomb angular integral - exchange part:

```math
b^{k}(lm_{l};l^{\prime}m_{l^{\prime}})=(2l+1)(2l^{\prime}+1)
\left(\begin{array}{ccc}
l & k & l^{\prime}\\
0 & 0 & 0
\end{array}\right)^{2}\left(\begin{array}{ccc}
l & k & l^{\prime}\\
-m_{l} & (m_{l}-m_{l^{\prime}}) & m_{l^{\prime}}
\end{array}\right)^{2}
```
#### Example:
```
b_exchange(1,1,1,2,2)
    2//5

b_exchange(6,3,2,3,-1)
    1050//20449
```
"""
function b_exchange(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)

    Base.iseven(k+l+l′) || return 0
    abs(l-l′) ≤ k ≤ l+l′ || return 0

    Δ = CamiMath.triangle_coefficient(l, k, l′)
    T = CamiMath._Racah_sqrt2(l, 0, k, 0, l′, 0) * CamiMath._Racah_sqrt2(l, -ml, k, ml-ml′, l′, ml′)
    S = CamiMath._Racah_sum(l, 0, k, 0, l′) * CamiMath._Racah_sum(l, -ml, k, ml-ml′, l′)

    o = abs((2l+1) * (2l′+1) * Δ * Δ * S * S * T)
    num = Int(numerator(o))
    den = Int(denominator(o))

    o = num//den

    return o

end

# ------------------------------------------------------------------------------
#           A_direct(k::Int, l::Int)
# ------------------------------------------------------------------------------

@doc raw"""
    A_direct(k::Int, l::Int)

Shell-averaged Coulomb angular integral - direct part:

```math
A^{k}(l)=\frac{1}{w_{l}}\sum_{m_{l}=-l}^{l}2a^{k}(lm_{l};l^{\prime}m_{l^{\prime}})=\frac{2l+1}{4l+1}\left(\begin{array}{ccc}
l & k & l\\
0 & 0 & 0
\end{array}\right)^{2}\!\!=\delta_{k,0},
```
where  ``w_{l}=2(2l+1)`` is the number of electrons in the closed shell ``(nl)^{w_l}``. Note that the shell average ``A^{k}(l)`` is *independent of* ``l``.
#### Example:
```
julia> l=2; wl=2(2l+1);

julia> [A_direct(k,l) for k=0:2l] == [1,0,0,0,0]
true

julia> [sum([2a_direct(k, l, ml, l, 0) for ml=-l:l])//wl for k=0:4] == [1,0,0,0,0]
true
```
"""
function A_direct(k::Int, l::Int)

    o = iszero(k) ? 1 : 0

    return o

end

# ------------------------------------------------------------------------------
#           B_exchange(k::Int, l::Int, l′::Int)
# ------------------------------------------------------------------------------

@doc raw"""
    B_exchange(k::Int, l::Int, l′::Int)

Shell-averaged Coulomb angular integral - exchange part:

```math
B^{k}(l,l^{\prime})=\frac{1}{w_{l}}\sum_{m_{l}=-l}^{l}b^{k}(lm_{l};l^{\prime}m_{l^{\prime}})=\tfrac{1}{2}\left(\begin{array}{ccc}
l & k & l^{\prime}\\
0 & 0 & 0
\end{array}\right)^{2}
```
where  ``w_{l}=2(2l+1)`` is the number of electrons in the closed shell ``(nl)^{w_l}``. Note that the shell average ``B^{k}(l,l^{\prime})`` is *independent of* ``m_{l^{\prime}}``.
#### Example:
```
julia> l = 2; l′= 3; wl = 2(2l+1);

julia> [B_exchange(k, l, l′) for k=0:(l+l′)] == [0, 3//70, 0, 2//105, 0, 5//231]
true

julia> [sum([b_exchange(k, l, ml, l′, 0) for ml=-l:l])//wl for k=0:(l+l′)] == [0, 3//70, 0, 2//105, 0, 5//231]
true
```
"""
function B_exchange(k::Int, l::Int, l′::Int)

    o = rationalize(CamiMath.threeJsymbol(l,0, k, 0, l′, 0; msg=false)^2)//2

    return o

end



# ======================== UGk(k, Z, grid) ===================================

@doc raw"""
    UGk(k::Int, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real

``k^{th}``-order-multipole contribution to the *exchange* screening potential 
of the (reduced) electronic wavefunctions `P1` and `P2` of the same atom.

```math
U_{G}^{k}(\rho)
=\frac{1}{\rho^{k+1}}\int_{0}^{\rho}\varrho^{k}\tilde{R}_{nl}(\varrho)
\tilde{R}_{n^{\prime}l^{\prime}}(\varrho)\,
\varrho^{2}d\varrho+\rho^{k}\int_{\rho}^{\infty}
\frac{1}{\varrho^{k+1}}\tilde{R}_{nl}(\varrho)
\tilde{R}_{n^{\prime}l^{\prime}}(\varrho)\,\varrho^{2}d\varrho.
```
#### Example:
```
atom = castAtom(Z=2, A=4, Q=0; msg=true)
orbit1 = castOrbit(n=1, ℓ=0; msg=true)
orbit2 = castOrbit(n=2, ℓ=0; msg=true)
scr = nothing
grid = autoGrid(atom, [orbit1,orbit2], Float64; Nboost=1, msg=true)
def1 = castDef(grid, atom, spinorbit1, codata; scr)
E = inEE(def1)
adams = castAdams(E, grid, def1)
E, def, adams, Z1 = adams_moulton_master(E, grid, def1, adams; Δν=Value(1,"kHz"), imax=50, msg=false);

def2 = castDef(grid, atom, spinorbit2, codata; scr)
E = inEE(def2)
adams = castAdams(E, grid, def2)
E, def, adams, Z2 = adams_moulton_master(E, grid, def2, adams; Δν=Value(1,"kHz"), imax=50, msg=false);

P1 = real(Z1);
P2 = real(Z2);

UG0 = UG(0, P1, P2, grid);
plot_function(UG0, 1:grid.N, grid; title="He4(1s,2s):  exchange screening potential")
```
The plot is made using `CairomMakie`.
NB.: `plot_function` is not included in the `CamiXon` package.
![Image](../assets/He41s-UG0.png)
"""
function UGk(k::Int, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real

    N = grid.N
    r = grid.r
    one = T(1)

    P1xP2 = P1 .* P2 
       
    inner = r.^k .* P1xP2  
    inner = [CamiDiff.grid_integration(inner, grid, 1:n) for n=1:N]
    UGk_inner = inner .* r[1:N].^-(k+1)
    UGk_inner = CamiDiff.regularize!(UGk_inner; k=3)

    outer = (one ./ r).^(k+1) .* P1xP2
    outer = CamiDiff.regularize!(outer; k=3)
    outer = [CamiDiff.grid_integration(outer, grid, n:N) for n=1:N]
    UGk_outer = outer .* r[1:N].^k 
    
    UGk = UGk_inner .+ UGk_outer

    return UGk

end

# ======================== UFk(k, Z, grid) ===================================

@doc raw"""
    UFk(k::Int, P::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real

``k^{th}``-order-multipole contribution to the *direct* screening potential 
by an electron in the (reduced) radial wavefunction `P` of an atom.

```math
U_{F}^{k}(\rho)
=\frac{1}{\rho^{k+1}}\int_{0}^{\rho}\varrho^{k}
\left[\tilde{R}_{nl}(\varrho)\right]^{2}
\varrho^{2}d\varrho+\rho^{k}\int_{\rho}^{\infty}
\frac{1}{\varrho^{k+1}}
\left[\tilde{R}_{nl}(\varrho)\right]^{2}\varrho^{2}d\varrho.
```
#### Example:
```
codata = castCodata(2022)
atom = castAtom(Z=2, A=4, Q=0; msg=false);
orbit = castOrbit(n=1, ℓ=0; msg=false);
grid = autoGrid(atom, orbit, Float64; msg=true);
def = castDef(grid, atom, orbit, codata);
E = 0;
scr = zeros(grid.T,grid.N);       
def, adams, inE, Z = adams_moulton_nodes(E, scr, grid, def; imax=100, msg=false);
def, adams, inE, Z = adams_moulton_iterate!(Z, inE, grid, def, adams; imax=25, ϵ=1e-10, msg=false);
P1 = real(Z);
UF0P1 = UF(0, P1, grid);
plot_function(scrUF0P1, 1:grid.N, grid; title="He4(1s,1s):  direct screening potential")
```
The plot is made using `CairomMakie`.
NB.: `plot_function` is not included in the `CamiXon` package.
![Image](../assets/He41s-UF0.png)
"""
function UFk(k::Int, P::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real

    return UGk(k, P, P, grid)

end

# ========================  UF(spinorbit1, spinorbit2, P, grid) ===================================

@doc raw"""
    UF(spinorbit1::Spinorbit, spinorbit2::Spinorbit, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real

Potential of *direct* screening for the spectator electron in `orbit 1` by the screening 
electron in `orbit2` with (reduced) radial wavefunction `P2`.

```math
U_{F}(u_{\kappa},u_{\kappa^{\prime}};\rho)
={\textstyle \sum\limits_{k=0}^{\infty}}a^{k}(lm_{l};l^{\prime}m_{l^{\prime}})U_{F}^{k}(nl;\rho).
```
#### Example:
```
```
"""
function UF(spinorbit1::Spinorbit, spinorbit2::Spinorbit, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
    
    l = spinorbit1.ℓ
    ml= spinorbit1.mℓ    
    l′ = spinorbit2.ℓ
    ml′= spinorbit2.mℓ
    kmax = 2 * min(l,l′)
    
    a = [a_direct(k, l, ml, l′, ml′) for k=0:2:kmax]
    potUF = [UFk(k, P2, grid) for k=0:2:kmax]
    
    UF = sum([a[i] .* potUF[i] for i ∈ eachindex(a)])
    UF = convert(Vector{T}, UF)
    
    return UF
    
end

# ========================  UG(spinorbit1, spinorbit2, P1, P2, grid) ===================================

@doc raw"""
    UG(spinorbit1::Spinorbit, spinorbit2::Spinorbit, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real

Potential of *exchange* screening of two electrons with (reduced) wavefunctions `P1` and `P2`,
corresponding to the electronic orbitals `orbit1` and `orbit2`.

NB.  `UG(spinorbit1, spinorbit2, P1, P2, grid)` = `UG(spinorbit2, spinorbit1, P1, P2, grid)`
```math
U_{G}(u_{\kappa},u_{\kappa^{\prime}};\rho)
={\textstyle \sum\limits_{k=0}^{\infty}}b^{k}(lm_{l};l^{\prime}m_{l^{\prime}})U_{G}^{k}(nl,n^{\prime}l^{\prime};\rho).
```
#### Example:
```
```
"""
function UG(spinorbit1::Spinorbit, spinorbit2::Spinorbit, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
    
    l = spinorbit1.ℓ
    ml= spinorbit1.mℓ    
    l′ = spinorbit2.ℓ
    ml′= spinorbit2.mℓ
    #kmax = 2 * min(l,l′)
    kmin = abs(l-l′)
    kmax = l+l′
    
    b = [b_exchange(k, l, ml, l′, ml′) for k=kmin:2:kmax]
    potUG = [UGk(k, P1, P2, grid) for k=0:2:kmax]
    
    UG = sum([b[i] .* potUG[i] for i ∈ eachindex(b)])
    UG = convert(Vector{T}, UG)
    
    return UG
    
end

# ------------------------------------------------------------------------------
#                       Fk(k, P1, P2, grid)
# ------------------------------------------------------------------------------

@doc raw"""   
    Fk(k::Int, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid) where T<:Real

``k^{th}``-order-multipole contribution to the *direct* radial integral over 
the (reduced) radial wavefunctions `P1` and `P2` of two electrons (in the orbitals 
``nl`` and ``n^{\prime}l^{\prime}``) in a central potential field.

```math
F^{k}(nl;n^{\prime}l^{\prime})
=\int_{0}^{\infty}U_{F}^{k}(nl;\rho)\left[\tilde{R}_{n^{\prime}l^{\prime}}(\rho)\right]^{2}\rho^{2}d\rho
=\int_{0}^{\infty}U_{F}^{k}(n^{\prime}l^{\prime};\rho)\left[\tilde{R}_{nl}(\rho)\right]^{2}\rho^{2}d\rho.
```

    Fk(k::Int, P::Vector{T}, grid::CamiDiff.Grid) where T<:Real

``k^{th}``-order contribution to the *direct* radial integral over the (reduced) 
radial wavefunction `P` of two *equivalent* ``nl`` electrons in a central potential.

```math
F^{k}(nl)
=\int_{0}^{\infty}U_{F}^{k}(nl;\rho)\left[\tilde{R}_{nl}(\rho)\right]^{2}\rho^{2}d\rho.
```
"""
function Fk(k::Int, P::Vector{T}, grid::CamiDiff.Grid) where T<:Real

    potUF = UFk(k, P, grid)
    PxP = P .* P

    Fk = CamiDiff.grid_integration(potUF .* PxP, grid)

    return Fk
    
end
function Fk(k::Int, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid) where T<:Real

    UF1 = UFk(k, P1, grid)
    P2xP2 = P2 .* P2

    Fk = CamiDiff.grid_integration(UF1 .* P2xP2, grid)

    return Fk
    
end

# ------------------------------------------------------------------------------
#                       Gk(k, P1, P2, grid)
# ------------------------------------------------------------------------------

@doc raw"""
    Gk(k::Int, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid) where T<:Real

``k^{th}``-order-multipole contribution to the *exchange* radial integral over 
the (reduced) radial wavefunctions `P1` and 'P2' of two electrons in a central potential.

```math
F^{k}(nl;n^{\prime}l^{\prime})
=\int_{0}^{\infty}U_{F}^{k}(nl;\rho)\left[\tilde{R}_{n^{\prime}l^{\prime}}(\rho)\right]^{2}\rho^{2}d\rho.
```
"""
function Gk(k::Int, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid) where T<:Real

    potUG = UGk(k, P1, P2, grid)
    P1xP2 = P1 .* P2

    Gk = CamiDiff.grid_integration(potUG .* P1xP2, grid)

    return Gk
    
end

@doc raw"""
    𝒥(spinorbit1::Spinorbit, spinorbit2::Spinorbit, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real

The *direct* integral of the electrostatic repulsion energy between two electrons in the (reduced) eigenstates `P1` and `P2` of an atom, 
which correspond to the orbitals `orbit1` and `orbit2`.

```math
\mathcal{J}(u_{\kappa},u_{\kappa^{\prime}})=(u_{\kappa},u_{\kappa^{\prime}}|\frac{1}{\rho_{12}}|u_{\kappa},u_{\kappa^{\prime}})
=\int_{0}^{\infty}U_{F}(u_{\kappa},u_{\kappa^{\prime}};\rho)\left[\tilde{R}_{n^{\prime}l^{\prime}}(\rho)\right]^{2}\!\rho^{2}d\rho.
```
"""
function 𝒥(spinorbit1::Spinorbit, spinorbit2::Spinorbit, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
    
    l = spinorbit1.ℓ
    ml= spinorbit1.mℓ    
    l′ = spinorbit2.ℓ
    ml′= spinorbit2.mℓ
    kmax = 2 * min(l,l′)
    
    potUF = UF(spinorbit1, spinorbit2, P2, grid)
    P2xP2 = P2 .* P2

    𝒥 = CamiDiff.grid_integration(potUF .* P2xP2, grid)
    
    return 𝒥
    
end

@doc raw"""
    𝒦(spinorbit1::Spinorbit, spinorbit2::Spinorbit, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real

The *exchange* integral

```math
\mathcal{K}(u_{\kappa},u_{\kappa^{\prime}})=(u_{\kappa},u_{\kappa^{\prime}}|\frac{1}{\rho_{12}}|u_{\kappa^{\prime}},u_{\kappa})
=\int_{0}^{\infty}U_{G}(u_{\kappa},u_{\kappa^{\prime}};\rho)\tilde{R}_{nl}(\rho)\tilde{R}_{n^{\prime}l^{\prime}}(\rho)\rho^{2}d\rho.
```
"""
function 𝒦(spinorbit1::Spinorbit, spinorbit2::Spinorbit, P1::Vector{T}, P2::Vector{T}, grid::CamiDiff.Grid{T}) where T<:Real
    
    l = spinorbit1.ℓ
    ml= spinorbit1.mℓ    
    l′ = spinorbit2.ℓ
    ml′= spinorbit2.mℓ
    
    kmin = abs(l-l′)
    kmax = l+l′
    
    potUG = UG(spinorbit1, spinorbit2, P1, P2, grid)
    P1xP2 = P1 .* P2

    𝒦 = CamiDiff.grid_integration(potUG .* P1xP2, grid)
    
    return 𝒦
    
end
