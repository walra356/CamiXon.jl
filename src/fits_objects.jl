# ...................................................... FITS objects .........................................................

"""
    FITS_name

FITS object to decompose the names of .fits files.

The fields are:
* `.name::String`: for filename 'p#.fits' this is 'p#.fits'
* `.prefix::String`: for filename 'p#.fits' this is 'p'
* `.numerator::String`: for filename 'p#.fits' this is '#', a serial number (e.g., '3') or a range (e.g., '3-7')
* `.extension::String`:  for filename 'p#.fits' this is '.fits'.
"""
struct FITS_name

    name::String
    prefix::String
    numerator::String
    extension::String

end

"""
    FITS_HDU

Object to hold a single "Header-Data Unit" (HDU).

The fields are
* `.filename::String`: name of the corresponding FITS file
* `.hduindex::Int`: identifier (a file may contain more than one HDU)
* `.header::T, where T=FITS_header`: the header object
* `.dataobject::V, where V=FITS_data`: the data object
"""
struct FITS_HDU{T,V}

    filename::String
    hduindex::Int
    header::T
    dataobject::V

end

"""
    FITS_header

Object to hold the header information of a `FITS_HDU`.

The fields are:
* `.hduindex::Int`: identifier (a file may contain more than one HDU)
* `.records::Array{String,1}`:  the header formated as an array of strings of 80 ASCII characters
* `.keys::Array{String,1}`: `keys[i]` - key corresponding to `records[i]` (record of index `i`)
* `.values::Array{Any,1}`: `value[i]` - corresponding to `records[i]`
* `.comments`: `comments[i]` - comment corresponding to `records[i]`
* `.dict::Dict{String,Any}`: Dictionary `key[i] => value[i]`
* `.maps::Dict{String,Int}`: Dictionary `key[i] => i`
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
    FITS_data

Object to hold the data of the `FITS_HDU` of given `hduindex` and `hdutype`.

The fields are:
* `.hduindex::Int`: identifier (a file may contain more than one HDU)
* `.hdutype::String`: accepted types are 'PRIMARY', 'IMAGE' and 'TABLE'
* `.data::Any`: in the from appropriate for the `hdutype`
"""
struct FITS_data

    hduindex::Int
    hdutype::String
    data

end

"""
    FITS_table

Object to hold the data of a `TABLE HDU` (a `FITS_HDU` for ASCII tables). It contains the data in the form of records (rows) of ASCII strings.

The fields are:
* `.hduindex::Int`: identifier (a file may contain more than one HDU)
* `.rows::Array{String,1}`: the table formated as an array of rows of ASCII strings
"""
struct FITS_table

    hduindex::Int
    rows::Array{String,1}

end

# ....................... parse FITS_TABLE into a Vector of its columns .........................................

"""
    parse_FITS_TABLE(hdu)

Parse `FITS_TABLE` into a Vector of its columns for further processing by the user.
#### Example:
```

f = fits_create("minimal.fits";protect=false)
fits_info(f[1])

 File: minimal.fits
 HDU: 1
 DataType: Any
 Datasize: (0,)

 Metainformation:
 SIMPLE  =                    T / file does conform to FITS standard
 NAXIS   =                    0 / number of data axes
 EXTEND  =                    T / FITS dataset may contain extensions
 COMMENT    Basic FITS file     / http://fits.gsfc.nasa.gov/iaufwg
 END

```
"""
function parse_FITS_TABLE(hdu::FITS_HDU)

    dict = hdu.header.dict
    thdu = Base.strip(Base.get(dict,"XTENSION", "UNKNOWN") ,['\'',' '])

    thdu == "TABLE" || return error("Error: $thdu is not an ASCII TABLE HDU")

    ncols = Base.get(dict,"TFIELDS", 0)
    nrows = Base.get(dict,"NAXIS2", 0)
    tbcol = [Base.get(dict,"TBCOL$n", 0) for n=1:ncols]
    tform = [Base.get(dict,"TFORM$n", 0) for n=1:ncols]
    ttype = [cast_FORTRAN_format(tform[n]).Type for n=1:ncols]
    tchar = [cast_FORTRAN_format(tform[n]).TypeChar for n=1:ncols]
    width = [cast_FORTRAN_format(tform[n]).width for n=1:ncols]
      itr = [(tbcol[k]:tbcol[k]+width[k]-1) for k=1:ncols]

     data = hdu.dataobject.data
     data = [[data[i][itr[k]] for i=1:nrows] for k=1:ncols]
     data = [tchar[k] == 'D' ? Base.join.(Base.replace!.(Base.collect.(data[k]), 'D'=>'E')) : data[k] for k=1:ncols]
     Type = [ttype[k] == "Aw" ? (width[k] == 1 ? Char : String) : ttype[k] == "Iw" ? Int : Float64 for k=1:ncols]
     data = [ttype[k] == "Aw" ? data[k] : parse.(Type[k],(data[k])) for k=1:ncols]

    return data

end
