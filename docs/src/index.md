```@meta
CurrentModule = CamiXon
```

# CamiXon

Here's an equation $\alpha$:

```math
\frac{n!}{k!(n - k)!} = \binom{n}{k}
```

This is the binomial coefficient.

---

To write a system of equations, use the `aligned` environment

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
