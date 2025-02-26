# Julia tools

```@docs
primitivetype(T::Type)
lc_primitivetype(o::Any)
lc_eltype(o)
conditionalType(n::T, nc::Int; msg=true) where {T<:Integer}
find_all(A::Union{String,AbstractArray{T,1}}, a::T...; count=false)  where T
find_first(A::Union{String,AbstractArray{T,1}}, a::T...; dict=false)  where T
find_last(A::Union{String,AbstractArray{T,1}}, a::T...; dict=false)  where T
```