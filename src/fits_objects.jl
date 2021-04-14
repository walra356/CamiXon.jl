# ...................................................... FITS objects .........................................................

struct FITS_name
    
    name::String
    prefix::String
    numerator::String
    extension::String
    
end

struct FITS_table
    
    hduindex::Int
    rows::Array{String,1}
    
end

struct FITS_header
    
    hduindex::Int
    records::Array{String,1}
    keys::Array{String,1}
    values::Array{Any,1}
    comments::Array{String,1}
    dict::Dict{String,Any}
    maps::Dict{String,Any}
    
end

struct FITS_data
    
    hduindex::Int
    hdutype::String
    data
    
end

struct FITS_HDU{T,V}
    
    filename::String    
    hduindex::Int
    header::T
    dataobject::V
    
end
