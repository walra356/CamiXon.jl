```@meta
CurrentModule = CamiXon
```

# CamiXon

Here's an equation:

```math
\frac{n!}{k!(n - k)!} = \binom{n}{k}
```

This is the binomial coefficient.

---

To write a system of equations, use the `aligned` environment:

```math
\begin{aligned}
\nabla\cdot\mathbf{E}  &= 4 \pi \rho \\
\nabla\cdot\mathbf{B}  &= 0 \\
\nabla\times\mathbf{E} &= - \frac{1}{c} \frac{\partial\mathbf{B}}{\partial t} \\
\nabla\times\mathbf{B} &= - \frac{1}{c} \left(4 \pi \mathbf{J} + \frac{\partial\mathbf{E}}{\partial t} \right)
\end{aligned}
```

These are Maxwell's equations.


```@docs
get_indices(A::AbstractArray{T,N}, a::T...)  where {T,N}
```

```@docs
get_indices_count(A::AbstractArray{T,N}, a::T...)  where {T,N}
```

```@docs
get_permutation_count(A::AbstractArray{T,N}; unique = false)  where {T,N}
```

```@index
```
