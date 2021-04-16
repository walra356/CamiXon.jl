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

# ........................................... cast filename into a FITSname object .................................

"""
    cast_FITS_name(filename::String)

Decompose the FITS filename 'filnam.fits' into its name, prefix, numerator and extension.
#### Examples:
```
strExample = "T23.01.fits"
f = analyze_FITS_name(strExample)
FITS_name("T23.01", "T23.", "01", ".fits")

f.name, f.prefix, f.numerator, f.extension
("T23.01", "T23.", "01", ".fits")
```
"""
function cast_FITS_name(str::String)

    Base.length(Base.strip(str)) == 0 && return error("Error: filename required")
    
    ne = Base.findlast('.',str)                                     # ne: first digit of extension
    nl = Base.length(str)                                           # ne: length of file name including extension
 
    hasextension = ne == nothing ? false : true

    if hasextension
        strNam = str[1:ne-1]
        strExt = Base.rstrip(str[ne:nl])
        strExt = Base.Unicode.lowercase(strExt)  
        isfits = strExt == ".fits" ? true : false
        n = Base.Unicode.isdigit(str[ne-1]) ? ne-1 : nothing        # n: last digit of numerator (if existent)
    else
        isfits = false
        n = Base.Unicode.isdigit(str[nl]) ? nl : nothing            # n: last digit of numerator (if existent)
    end
    
    isfits || return println("FitsError: '$str': filename lacks mandatory '.fits' extension")

    if n != nothing
        strNum = ""
        while Base.Unicode.isdigit(str[n])
            strNum = str[n] * strNum
            n -= 1
        end
        strPre = str[1:n]
    else
        strPre = strNam
        strNum = " "
    end

    return FITS_name(strNam,strPre,strNum,strExt)
    
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
* `.keys::Array{String,1}`: keys[i] - key corresponding to `records[i]` for record index `i` 
* `.values::Array{Any,1}`: value[i] - corresponding to `records[i]`
* `.comments`: comments[i] - comment corresponding to given `records[i]`
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

# ........................................... cast records into a FITS_header object .................................

function _cast_header(records::Array{String,1}, hduindex::Int)
   
    records = _rm_blanks(records)         # remove blank records to collect header records data (key, val, comment) 
    nrec = length(records)                # number of keys in header with given hduindex

    keys = [Base.strip(records[i][1:8]) for i=1:nrec]
    vals = [records[i][9:10] ≠ "= " ? records[i][11:31] : _fits_parse(records[i][11:31]) for i=1:nrec]
    coms = [records[i][34:80] for i=1:nrec]
    dict = [keys[i] => vals[i] for i=1:nrec]
    maps = [keys[i] => i for i=1:nrec]
    
    return FITS_header(hduindex,records,keys,vals,coms,Dict(dict),Dict(maps))
    
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

# .................................. cast data into FITS_data objects ....................................

function _cast_data(hduindex::Int, hdutype::String, data)
    
    return FITS_data(hduindex, hdutype, data)

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
    parse_FITS_TABLE(FITS_HDU)

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
function parse_FITS_TABLE(FITS_HDU)
    
    dict = FITS_HDU.header.dict
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
    
     data = FITS_HDU.dataobject.data
     data = [[data[i][itr[k]] for i=1:nrows] for k=1:ncols]
     data = [tchar[k] == 'D' ? Base.join.(Base.replace!.(Base.collect.(data[k]), 'D'=>'E')) : data[k] for k=1:ncols]
     Type = [ttype[k] == "Aw" ? (width[k] == 1 ? Char : String) : ttype[k] == "Iw" ? Int : Float64 for k=1:ncols]
     data = [ttype[k] == "Aw" ? data[k] : parse.(Type[k],(data[k])) for k=1:ncols]
    
    return data
    
end
