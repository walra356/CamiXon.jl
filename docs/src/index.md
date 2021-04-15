```@meta
CurrentModule = CamiXon
```

# CamiXon.jl

A package for image analysis of backscattered light

---
## Table of contents

```@contents
```

## FITS

```@docs
FITS_HDU
FITS_table
parse_FITS_TABLE(FITS_HDU)
FITS_header
FITS_data
FITS_name
cast_FITS_name(filename::String)
print_hduinfo(FITS_HDU)
fits_create(filename::String, data=[]; protect=true)
fits_read(filename::String)
fits_extend(filename::String, data_extend, hdutype="IMAGE") 
fits_copy(filenameA::String, filenameB::String=" "; protect=true)
fits_add_key(filename::String, hduindex::Int, key::String, val::Real, com::String)
fits_edit_key(filename::String, hduindex::Int, key::String, val::Real, com::String)
fits_delete_key(filename::String, hduindex::Int, key::String)
```

## FORTRAN 

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
