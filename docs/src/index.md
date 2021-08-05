```@meta
CurrentModule = CamiXon
```

# CamiXon.jl

A package for image analysis of backscattered light

---
## Table of contents

```@contents
```
## Finite-difference methods

Consider an analytic function ``f`` tabulated in *standard ordering of growing index* at ``n`` positions on a *grid*. The *finite difference* of two adjacent values on a *uniform grid* is given by the relation

```math
\nabla f[n] = f[n]-f[n-1].
```

This is called the *backward difference* notation. In this notation the  ``k^{th}``-*order finite differences* (``k+1``-point finite differences) are defined given by a *weighted sum* over the function values ``f[n],\ \ldots,\ f[n-k]``,

```math
\nabla^k f[n] = f[n] + c_1^kf[n-1] + \cdots + c_k^kf[n-k] = \sum_{j=0}^{k} c_j^kf[n-j] = \sum_{j=0}^{k} c_{k-j}^kf[n-k+j].
```

The k+1 coefficients ``c_{j}^{k}=(-1)^{j}\binom{k}{j}`` are *weight factors* (short: *weights*) defining the summation. Note that ``c_{0}^{k}\equiv1$ and $c_{k}^{k}=(-1)^{k}$``.

Turning to the standard ordering of terms the summation becomes

```math
\nabla^k f[n] = \sum_{j=0}^{k} c_{k-j}^kf[n-k+j].

Functions:  

`f_diff_weight(k,i)` `` \rightarrow c_i^k``

`f_diff_weights(k)` `` \rightarrow \ [c_k^k,\ c_1^k,\ldots,\ c_0^k]$``

`f_diff_weights_array(kmax)` `` \rightarrow \ [\ [c_0^0],\ [c_1^1,c_0^1],\ \ldots,\ [c_k^k,\ c_{k-1}^k,\ldots,\ c_0^k] ]``


```@docs
f_diff_weight(k::Int, i::Int)
f_diff_weights(k::Int)
f_diff_weights_array(kmax::Int)
f_diff_expansion_weights(coeffs, ∇)
f_diff_expansion_coeffs_interpolation(k::Int, x::T) where T<:Real
interpolation_offset_positions(n::Int, k::Int, i::Int)
summation_range(n::Int, j::Int, k::Int, i::Int)
```

## FITS

FITS stands for 'Flexible Image Transport System'. This is an open standard origionally developed for the astronomy community to store telescope images together with tables of spectral information. Over the years it has developed into a scientific standard - http://fits.gsfc.nasa.gov/iaufwg.

Within CamiXion only the basic FITS functionality is implemented for users not requiring celestal coordinates. The user can create, read and extend .fits files as well as create, edit and delete user-defined metainformation.

A FITS file consists of a sequence of one or more header-data-units (HDUs), each containing a data block preceeded by header records of metainformation.

By the command `f = fits_read(filnam)` we asign a collection of `FITS_HDU` objects from the file `filnam` to the variable `f`.

### FITS - Types

```@docs
FITS_HDU
FITS_header
FITS_data
FITS_table
FITS_name
```

### FITS - HDU Methods

```@docs
fits_info(hdu::FITS_HDU)
parse_FITS_TABLE(hdu::FITS_HDU)
```

### FITS - File Methods

```@docs
cast_FITS_name(filename::String)
fits_combine(filnamFirst::String, filnamLast::String; protect=true)
fits_copy(filenameA::String, filenameB::String=" "; protect=true)
fits_create(filename::String, data=[]; protect=true)
fits_extend(filename::String, data_extend, hdutype="IMAGE")
fits_read(filename::String)
```

### FITS - Key Methods

```@docs
fits_add_key(filename::String, hduindex::Int, key::String, val::Real, com::String)
fits_delete_key(filename::String, hduindex::Int, key::String)
fits_edit_key(filename::String, hduindex::Int, key::String, val::Real, com::String)
fits_rename_key(filename::String, hduindex::Int, keyold::String, keynew::String)
```

## FORTRAN

```@docs
FORTRAN_format
cast_FORTRAN_format(str::String)
cast_FORTRAN_datatype(str::String)
```

## Plotting

```@docs
step125(x::Real)
select125(x)
steps(x::Vector{T} where T<:Real)
stepcenters(x::Vector{T} where T<:Real)
stepedges(x::Vector{T} where T<:Real)
edges(px, Δx=1.0, x0=0.0)
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
log10_characteristic_power(x)
log10_mantissa(x)
polynom_deriv_coeffs(c,deriv=0)
polynom(c::Vector{T}, x::T) where T<:Real
```

## Index

```@index
```
