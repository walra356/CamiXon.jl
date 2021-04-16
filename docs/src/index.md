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

FITS stands for 'Flexible Image Transport System'. This is an open standard origionally developed for the astronomy community to store telescope images together with tables of spectral information. Over the years it has developed into a scientific standard - http://fits.gsfc.nasa.gov/iaufwg.

Within CamiXion only the basic FITS functionality is implemented for users not requiring celestal coordinates. The user can create, read and extend .fits files as well as create, edit and delete user-defined metainformation.

A FITS file consists of a sequence of one or more header-data-units (HDUs), each containing a data block preceeded by header records of metainformation.

By the command `f = fits_read(filnam)` we asign a collection of `FITS_HDU` objects from the file `filnam` to the variable `f`.


```@docs
FITS_HDU
FITS_header
FITS_data
FITS_table
parse_FITS_TABLE(FITS_HDU)
FITS_name
cast_FITS_name(filename::String)
fits_info(FITS_HDU)
test_fits_info(FITS_HDU)
fits_create(filename::String, data=[]; protect=true)
test_fits_create()
fits_read(filename::String)
test_fits_read()
fits_extend(filename::String, data_extend, hdutype="IMAGE") 
test_fits_extend() 
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
