# ...................................................... FITS objects .........................................................

"""
    FITS_HDU

Object to hold a single "Header-Data Unit" (HDU) of a give FITS file. 

Object to hold `FITS_header` and `FITS_data` objects of a `FITS_HDU` of given `hduindex` and `hdutypes` from file 'filename'. 

The fields are 
* `.filename::String`
* `.hduindex::Int`: identifier (a file may contain more than one HDU)
* `.header::T, where T=FITS_header`
* `.dataobject::V, where V=FITS_data` 
"""
struct FITS_HDU{T,V}
    
    filename::String    
    hduindex::Int
    header::T
    dataobject::V
    
end


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
    FITS_table

Object to hold the rows of an ASCII `TABLE HDU` of given 'hduindex'. 

The fields are:
* `.hduindex::Int`: identifier (a file may contain more than one HDU)
* `.rows::Array{String,1}`: the table formated as an array of rows os ASCII strings
"""
struct FITS_table
    
    hduindex::Int
    rows::Array{String,1}
    
end

"""
    FITS_header

Object to hold the header information of a `FITS_HDU` of given `hduindex`. 

The fields are:
* `.hduindex::Int`: identifier (a file may contain more than one HDU)
* `.records::Array{String,1}`:  the header formated as an array of strings of 80 ASCII characters
* `.keys::Array{String,1}`: keys[i] - key corresponding to `records[i]`
* `.values::Array{Any,1}`: value[i] - corresponding to `records[i]
* `.comments`: comments[i] - comment corresponding to given `records[i]
* `.dict::Dict{String,Any}`: Dictionary key[i] => value[i]
* `.maps::Dict{String,Int}`: Dictionary key[i] => i
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

# .................................. cast into FITS_header and FITS_data objects ....................................

function _cast_data(hduindex::Int, hdutype::String, data)
    
    return FITS_data(hduindex, hdutype, data)

end

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

function _rm_blanks(records::Array{String,1})               # remove blank records
    
    record_B = repeat(' ',length(records[1]))
    
    return [records[i] for i ∈ findall(records .≠ record_B)]
    
end
