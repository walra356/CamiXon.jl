# ...................................................... FITS_HDU Objects .........................................................

"""
    FITS_HDU{T,V}

Object to hold a single "Header-Data Unit" (HDU).

The fields are
* `.filename::String`:  name of the corresponding FITS file
* `   .hduindex::Int`:  identifier (a file may contain more than one HDU)
* `       .header::T`:  the header object where T=FITS_header
* `   .dataobject::V`:  the data object where V=FITS_data
"""
struct FITS_HDU{T,V}

    filename::String
    hduindex::Int
    header::T        #FITS_header
    dataobject::V    #FITS_data

end

# ........................................... FITS_name Object..........................................................

"""
    FITS_name

FITS object to decompose the names of .fits files.

The fields are:
* `     .name::String`:  for 'p#.fits' this is 'p#.fits'
* `   .prefix::String`:  for 'p#.fits' this is 'p'
* `.numerator::String`:  for 'p#.fits' this is '#', a serial number (e.g., '3') or a range (e.g., '3-7')
* `.extension::String`:  for 'p#.fits' this is '.fits'.
"""
struct FITS_name

    name::String
    prefix::String
    numerator::String
    extension::String

end

# ........................................... FITS_header Object..........................................................

"""
    FITS_header

Object to hold the header information of a `FITS_HDU`.

The fields are:
* `           .hduindex::Int`:  identifier (a file may contain more than one HDU)
* `.records::Array{String,1}`:  the header formated as an array of strings of 80 ASCII characters
* `   .keys::Array{String,1}`:  `keys[i]` - key corresponding to `records[i]` (record of index `i`)
* `    .values::Array{Any,1}`:  `value[i]` - corresponding to `records[i]`
* `        .comments::String`:  `comments[i]` - comment corresponding to `records[i]`
* `  .dict::Dict{String,Any}`:  Dictionary `key[i] => value[i]`
* `  .maps::Dict{String,Int}`:  Dictionary `key[i] => i`
"""
struct FITS_header

    hduindex::Int
    records::Array{String,1}
    keys::Array{String,1}
    values::Array{Any,1}
    comments::Array{String,1}
    dict::Dict{String,Any}
    maps::Dict{String,Int}

end

# ........................................... FITS_data Object ...................................................

"""
    FITS_data

Object to hold the data of the `FITS_HDU` of given `hduindex` and `hdutype`.

The fields are:
* `  .hduindex::Int`:  identifier (a file may contain more than one HDU)
* `.hdutype::String`:  accepted types are 'PRIMARY', 'IMAGE' and 'TABLE'
* `      .data::Any`:  in the from appropriate for the `hdutype`
"""
struct FITS_data

    hduindex::Int
    hdutype::String
    data

end

# ........................................... FITS_table Object ..........................................................

"""
    FITS_table

Object to hold the data of a `TABLE HDU` (a `FITS_HDU` for ASCII tables). It contains the data in the form of records (rows) of ASCII strings.

The fields are:
* `        .hduindex::Int`:  identifier (a file may contain more than one HDU)
* `.rows::Array{String,1}`:  the table formated as an array of rows of ASCII strings
"""
struct FITS_table

    hduindex::Int
    rows::Array{String,1}

end
