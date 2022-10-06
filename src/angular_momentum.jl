# ========= threeJsymbol(j1, m1, j2, m2, j3, m3; msg=false) ====================

# ........................................................
function _threeJroot2(j1, m1, j2, m2, J, M)

    a = Int(j1 + m1)
    b = Int(j1 - m1)
    c = Int(j2 + m2)
    d = Int(j2 - m2)
    e = Int(J + M)
    f = Int(J - M)

    a = a ≥ 0 ? factorialbig(a) : return 0
    b = b ≥ 0 ? factorialbig(b) : return 0
    c = c ≥ 0 ? factorialbig(c) : return 0
    d = d ≥ 0 ? factorialbig(d) : return 0
    e = e ≥ 0 ? factorialbig(e) : return 0
    f = f ≥ 0 ? factorialbig(f) : return 0

    return a * b * c * d * e * f

end
# ........................................................
function _Racah_denom(j1, m1, j2, m2, J, t:: Int)

    a = Int(J - j2 + t + m1)
    b = Int(J - j1 + t - m2)
    c = Int(j1 + j2 - J - t)
    d = Int(j1 - t - m1)
    e = Int(j2 - t + m2)

    a = a ≥ 0 ? factorialbig(a) : return 0
    b = b ≥ 0 ? factorialbig(b) : return 0
    c = c ≥ 0 ? factorialbig(c) : return 0
    d = d ≥ 0 ? factorialbig(d) : return 0
    e = e ≥ 0 ? factorialbig(e) : return 0
    t = t ≥ 0 ? factorialbig(t) : return 0

    return a * b * c * d * e * t

end
# ........................................................
function _Racah_sum(j1, m1, j2, m2, J)

    o = 0
    t = 0

    while t < 25
        sign = iseven(t) ? 1 : -1
        d = _Racah_denom(j1, m1, j2, m2, J, t)
        d > 0 || break
        o += (sign//d)
        t += 1
    end

    return o

end
# ........................................................

@doc raw"""
    threeJsymbol(j1::Real, m1::Real, j2::Real, m2::Real, j3::Real, m3::Real; msg=false)

Wigner 3j symbol. This is a vector coupling coefficient with optimized symmetry
properties. The 3j symbols are zero unless ``Δ(j_{1},j_{2},j_{3})>0``
(triangle inequality holds) and ``m_{1}+m_{2}+m_{3}=0``. The implementation
is based on the Racah formula:

```math
\left(\begin{array}{ccc}
j_{1} & j_{2} & j_{3}\\
m_{1} & m_{2} & m_{3}
\end{array}\right)=
(-1)^{j_{1}-j_{2}-m_{3}}\sqrt{\Delta(j_{1}j_{2}J)}\\\times
\sqrt{\left(j_{1}+m_{1}\right)!
\left(j_{1}-m_{1}\right)!
\left(j_{2}+m_{2}\right)!
\left(j_{2}-m_{2}\right)!
\left(j_{3}+m_{3}\right)!
\left(j_{3}-m_{3}\right)!}
\\\times\sum_{t}\frac{(-)^{t}}{t!(j_{3}-j_{2}+t+m_{1})!
(j_{3}-j_{1}+t-m_{2})!
(j_{1}+j_{2}-j_{3}-t)!(j_{1}-t-m_{1})!(j_{2}-t+m_{2})!}
```
#### Example:
```
o = threeJsymbol(3, 0, 4, -1, 5, 1; msg=true); println(" = $o")
    -√(361/30030) = -0.10964174397241236

threeJsymbol(3, 0, 4, -1, 5, 1)
    -0.10964174397241236

threeJsymbol(0, 0, 0, 0, 0, 0)
    1.0
```
"""
function threeJsymbol(j1::Real, m1::Real, j2::Real, m2::Real, j3::Real, m3::Real; msg=false)

    (j1,m1,j2,m2,j3,m3) = promote(j1,m1,j2,m2,j3,m3)

    iszero(m1 + m2 + m3) || return 0

    Δ = triangle_coefficient(j1, j2, j3)
    T = _threeJroot2(j1, m1, j2, m2, j3, m3)
    R = _Racah_sum(j1, m1, j2, m2, j3)
    S = R * R
    A = Δ * T * S

    sgn_phase = iseven(j1-j2+m3) ? 1 : -1
    sgn_racah = sign(R)
    sgn = sgn_phase * sgn_racah

    msg && print((sgn < 0 ? "-" : "") * "√(" * strRational(A) * ")")

    return sgn * sqrt(A)

end

# =============== CGC(j1, m1, j2, m2l, J, M; msg=false) ========================

@doc raw"""
    CGC(j1::Real, m1::Real, j2::Real, m2::Real, J::Real, M::Real; msg=false)

Clebsch-Gordan coefficient (CGC). This is a vector-coupling coefficient in
Dirac notation. The CGCs are zero unless ``Δ(j_{1},j_{2},j_{3})>0``
(triangle inequality holds) and ``M=m_{1}+m_{2}``. The relation to the
Wigner 3j symbols is given by:

```math
\langle j_{1}m_{1};j_{2}m_{2}|JM\rangle\equiv
(-1)^{j_{1}-j_{2}+M}\sqrt{2J+1}\left(\begin{array}{ccc}
j_{1} & j_{2} & J\\
m_{1} & m_{2} & -M
\end{array}\right)
```
#### Example:
```
j1=3; m1=0
j2=4; m2=-1
J=5; M=-1
o = CGC(j1, m1, j2, m2, J, M; msg=true); println(" = $o")
o = CGC(j1, m1, j2, m2, J, M); println(o)
o = (-1)^(j1-j2+M) * sqrt(2J+1) * threeJsymbol(j1, m1, j2, m2, J, -M); println(o)
    -√(361/2730) = -0.36364052611670256
    -0.36364052611670256
    -0.36364052611670256
```
"""
function CGC(j1::Real, m1::Real, j2::Real, m2::Real, J::Real, M::Real; msg=false)

    (j1,m1,j2,m2,J,M) = promote(j1,m1,j2,m2,J,M)

    sgn = iseven(j1-j2+M) ? 1 : -1
    tJs = threeJsymbol(j1, m1, j2, m2, J, -M)

    if msg
        Δ = triangle_coefficient(j1, j2, J)
        T = _threeJroot2(j1, m1, j2, m2, J, M)
        R = _Racah_sum(j1, m1, j2, m2, J)
        S = R * R
        A = Δ * T * S
        s = sign(R) < 0 ? "-" : ""
        print(s * "√(" * strRational(A * (2J+1)) * ")")
    end

    return sgn * sqrt(2J+1) * tJs

end
