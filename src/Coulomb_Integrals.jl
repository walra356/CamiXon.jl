@doc raw"""
    a_coeff(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)

Angular coefficient for the *direct* Coulomb integral:

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
a(2,1,1,2,2)
    2//35

a(6,3,2,3,-1)
    -250//20449
```
"""
function a_coeff(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)

    Base.iseven(k) || return 0
    0 ≤ k ≤ 2min(l,l′) || return 0
    k > 0 || return 1

    Δ1 = triangle_coefficient(l, k, l)
    Δ2 = triangle_coefficient(l′, k, l′)
    T = _threeJroot2(l, 0, k, 0, l, 0) * _threeJroot2(l, -ml, k, 0, l, ml) * _threeJroot2(l′, 0, k, 0, l′, 0) * _threeJroot2(l′, -ml′, k, 0, l′, ml′)
    S = _Racah_sum(l, 0, k, 0, l) * _Racah_sum(l, -ml, k, 0, l) * _Racah_sum(l′, 0, k, 0, l′) * _Racah_sum(l′, -ml′, k, 0, l′)

    a = (2l+1) * (2l′+1) * Δ1 * Δ2 * S
    o = a * a * T
    num = Int(sqrt(numerator(o)))
    den = Int(sqrt(denominator(o)))

    o = sign(S) * num//den

    return o

end

@doc raw"""
    b_coeff(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)

Angular coefficient for the *exchange* Coulomb integral:

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
b(1,1,1,2,2)
    2//5

b(6,3,2,3,-1)
    1050//20449
```
"""
function b_coeff(k::Int, l::Int, ml::Int, l′::Int, ml′::Int)

    Base.iseven(k+l+l′) || return 0
    abs(l-l′) ≤ k ≤ l+l′ || return 0

    Δ = triangle_coefficient(l, k, l′)
    T = _threeJroot2(l, 0, k, 0, l′, 0) * _threeJroot2(l, -ml, k, ml-ml′, l′, ml′)
    S = _Racah_sum(l, 0, k, 0, l′) * _Racah_sum(l, -ml, k, ml-ml′, l′)

    o = abs((2l+1) * (2l′+1) * Δ * Δ * S * S * T)
    num = Int(numerator(o))
    den = Int(denominator(o))

    o = num//den

    return o

end
