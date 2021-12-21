# ...................................................... FITS objects .........................................................

"""
    FITS_name(name::String, prefix::String, numerator::String, extension::String)

FITS object to decompose the names of .fits files.

The fields are:
* `     .name`:  for 'p#.fits' this is 'p#.fits'
* `   .prefix`:  for 'p#.fits' this is 'p'
* `.numerator`:  for 'p#.fits' this is '#', a serial number (e.g., '3') or a range (e.g., '3-7')
* `.extension`:  for 'p#.fits' this is '.fits'.
"""
struct FITS_name

    name::String
    prefix::String
    numerator::String
    extension::String

end

"""
    FITS_HDU{T,V}(filename::String, hduindex::Int, header::T, dataobject::V) where T,V = FITS_header, FITS_data

Object to hold a single "Header-Data Unit" (HDU).

The fields are
* `  .filename`:  name of the corresponding FITS file
* `  .hduindex`:  identifier (a file may contain more than one HDU)
* `    .header`:  the header object
* `.dataobject`:  the data object
"""
struct FITS_HDU{T,V} where {T,V} = {FITS_header, FITS_data}

    filename::String
    hduindex::Int
    header::T
    dataobject::V

end

"""
    FITS_header(hduindex, records, keys, values, comments, dict, maps)

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

# ........................................... FITS_data object ..........................................................

"""
    FITS_data(hduindex::Int, hdutype::String, data::Any)

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

"""
    FITS_table(hduindex::Int, rows::Array{String,1})

Object to hold the data of a `TABLE HDU` (a `FITS_HDU` for ASCII tables). It contains the data in the form of records (rows) of ASCII strings.

The fields are:
* `.hduindex`:  identifier (a file may contain more than one HDU)
* `    .rows`:  the table formated as an array of rows of ASCII strings
"""
struct FITS_table

    hduindex::Int
    rows::Array{String,1}

end
