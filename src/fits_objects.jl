# ...................................................... FITS objects .........................................................

"""
    FITS_name

FITS object to decompose filenames of the type 'p#.fits'. The fields are '.name' (p#.fits), '.prefix' (p), '.numerator' (#) and '.extension' (.fits).
"""
struct FITS_name
    
    name::String
    prefix::String
    numerator::String
    extension::String
    
end

"""
    FITS_table

Object to hold the rows of an ASCII TABLE HDU of given 'hduindex'. The fields are '.hduindex' and '.rows'.
"""
struct FITS_table
    
    hduindex::Int
    rows::Array{String,1}
    
end

"""
    FITS_header

Object to hold the header information of a FITS_HDU of given 'hduindex'. 
The fields are '.hduindex', '.records', '.keys', '.values', '.comments' '.dict', '.maps'.
"""
struct FITS_header
    
    hduindex::Int
    records::Array{String,1}
    keys::Array{String,1}
    values::Array{Any,1}
    comments::Array{String,1}
    dict::Dict{String,Any}
    maps::Dict{String,Any}
    
end

"""
    FITS_data

Object to hold the data of a FITS_HDU of given 'hduindex' and 'hdutypes'. 
The accepted types are 'PRIMARY', 'IMAGE' and 'TABLE' ('BINTABLE' is not implemented). 
The fields are '.hduindex', '.hdutypes' and '.data'.
"""
struct FITS_data
    
    hduindex::Int
    hdutype::String
    data
    
end

"""
    FITS_HDU

Object to hold 'FITS_header' and 'FITS_data' objects of 'FITS_HDU' of given 'hduindex' and 'hdutypes' from file 'filename'. 
The fields are '.filename', '.hduindex', '.header' and '.dataobject'.
"""
struct FITS_HDU{T,V}
    
    filename::String    
    hduindex::Int
    header::T
    dataobject::V
    
end
