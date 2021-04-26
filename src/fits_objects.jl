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
```strExample = "example.fits"
data = [10, 20, 30]
fits_create(strExample, data; protect=false)

t1 = Float16[1.01E-6,2.0E-6,3.0E-6,4.0E-6,5.0E-6]
t2 = [0x0000043e, 0x0000040c, 0x0000041f, 0x0000042e, 0x0000042f]
t3 = [1.23,2.12,3.,4.,5.]
t4 = ['a','b','c','d','e']
t5 = ["a","bb","ccc","dddd","ABCeeaeeEEEEEEEEEEEE"]
data = [t1,t2,t3,t4,t5]
fits_extend(strExample, data, "TABLE")

f = fits_read(strExample)
d = f[2].header.dict
d = [get(d,"TFORM$i",0) for i=1:5]; println(strip.(d))
  SubString{String}["'E6.1    '", "'I4      '", "'F4.2    '", "'A1      '", "'A20     '"]

f[2].dataobject.data                            # this is the table hdu
  5-element Vector{String}:
   "1.0e-6 1086 1.23 a a                    "
   "2.0e-6 1036 2.12 b bb                   "
   "3.0e-6 1055 3.0  c ccc                  "
   "4.0e-6 1070 4.0  d dddd                 "
   "5.0e-6 1071 5.0  e ABCeeaeeEEEEEEEEEEEE "

parse_FITS_TABLE(f[2])
  5-element Vector{Vector{T} where T}:
   [1.0e-6, 2.0e-6, 3.0e-6, 4.0e-6, 5.0e-6]
   [1086, 1036, 1055, 1070, 1071]
   [1.23, 2.12, 3.0, 4.0, 5.0]
   ["a", "b", "c", "d", "e"]
   ["a                   ", "bb                  ", "ccc                 ", "dddd                ", "ABCeeaeeEEEEEEEEEEEE"]
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
