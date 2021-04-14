```@meta
CurrentModule = CamiXon
```

# CamiXon.jl

A package for image analysis of backscattered light

---
## Table of contents

```@contents
```

## File manipulation

```@docs
decompose_filnam(str::String)
fits_combine(filnamFirst::String, filnamLast::String; info=false)
fits_copy(filnam, filnamOut=""; protect=true)
fits_info(filnam::String; info=false)
fits_key_create(filnam::String, key::String, value, comment::String)
fits_key_delete(filnam::String, key::String)
fits_key_edit(filnam::String, key::String, value, comment::String)
fits_key_info(filnam::String, key::String)
fits_key_rename(filnam::String, key::String, keynew::String)
```
## FORTRAN format 

```@docs
FORTRAN_format
cast_FORTRAN_format(str::String)
cast_FORTRAN_datatype(str::String)
```

## Search algorithms

```@docs
find_all(A::Union{String,AbstractArray{T,1}}, a::T...; count=false)  where T
find_first(A::Union{String,AbstractArray{T,1}}, a::T...; dict=false)  where T
find_last(A::Union{String,AbstractArray{T,1}}, a::T...; dict=false)  where T
```

## Math

```@docs
canonical_partitions(n::Int, m=0; header=true, reverse=true)
integer_partitions(n::Int, m=0; transpose=false, count=false)
```

## Index

```@index
```
