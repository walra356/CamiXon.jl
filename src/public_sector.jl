# ......................................... FITS public sector .................................................................

"""
    fits_info(FITS_HDU)

Print metafinformation of given `FITS_HDU`
#### Example:
```

strExample = "minimal.fits"
f = fits_create(strExample; protect=false)
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
function fits_info(FITS_HDU)
    
    info = [
        "\r\nFile: " * FITS_HDU.filename,
        "hdu: " * string(FITS_HDU.hduindex),
        "hdutype: " * FITS_HDU.dataobject.hdutype,
        "DataType: " * string(Base.eltype(FITS_HDU.dataobject.data)),
        "Datasize: " * string(Base.size(FITS_HDU.dataobject.data)),
        "\r\nMetainformation:"
        ]
    
    records = FITS_HDU.header.records
    
    append!(info,_rm_blanks(records))

    return print(Base.join(info .* "\r\n"))
    
end

"""
    fits_create(filename [, data [; protect=true]])

Create FITS file of given filename [, optional data block [, default overwrite protection]] and return Array of HDUs
#### Examples:
```

f = fits_create("minimal.fits";protect=false)
a = f[1].dataobject.data
b = f[1].header.keys
println(a);println(b)
 Any[]
 ["SIMPLE", "NAXIS", "EXTEND", "COMMENT", "END"]

data = [0x0000043e, 0x0000040c, 0x0000041f]
f = fits_create("remove.fits", data; protect=false)
a = f[1].dataobject.data
b = f[1].header.keys
println(a);println(b)
 UInt32[0x0000043e, 0x0000040c, 0x0000041f]
 ["SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "BZERO", "BSCALE", "EXTEND", "COMMENT", "END"]

```
"""
function fits_create(filename::String, data=[]; protect=true)
    
    _validate_FITS_name(filename)                            
    _isavailable(filename, protect) || error("FitsError: '$filename': creation failed") 

    nhdu = 1
    hdutype = "PRIMARY"

       FITS_data = [_cast_data(i, hdutype, data) for i=1:nhdu] 
    FITS_headers = [_cast_header(_PRIMARY_input(FITS_data[i]), i) for i=1:nhdu]
    
    FITS = [FITS_HDU(filename, i, FITS_headers[i], FITS_data[i]) for i=1:nhdu]
    
    _fits_save(FITS)
       
    return FITS
    
end


"""
    fits_read(filename)

Read FITS file and return Array of `FITS_HDU`'s
#### Examples:
```

f = fits_create("minimal.fits";protect=false)
f[1].dataobject.data

 Any[]

g = fits_read("minimal.fits")
g[1].dataobject.data

 Any[]

```
"""
function fits_read(filename::String)
     
    o = _fits_read_IO(filename)
    
    nhdu = _hdu_count(o)
    
    FITS_headers = [_read_header(o,i) for i=1:nhdu] 
       FITS_data = [_read_data(o,i) for i=1:nhdu] 
    
    FITS = [FITS_HDU(filename, i, FITS_headers[i], FITS_data[i]) for i=1:nhdu]
        
    return FITS
    
end

"""
    fits_extend(filename, data_extend, hdutype="IMAGE") 

Extend the FITS file of given filename with the data of `hdutype` from `data_extend`  and return Array of HDUs
#### Examples:
```

f = fits_read("minimal.fits")

```
"""
function fits_extend(filename::String, data_extend, hdutype="IMAGE")   
        
    hdutype == "IMAGE"    ? (records, data) = _IMAGE_input(data_extend)    : 
    hdutype == "TABLE"    ? (records, data) = _TABLE_input(data_extend)    :  
    hdutype == "BINTABLE" ? (records, data) = _BINTABLE_input(data_extend) : error("FitsError: unknown HDU type")
    
    o = _fits_read_IO(filename)  
    
    nhdu = _hdu_count(o)
    
    FITS_headers = [_read_header(o,i) for i=1:nhdu] 
       FITS_data = [_read_data(o,i) for i=1:nhdu] 
    
    nhdu = nhdu + 1
    
    push!(FITS_headers, _cast_header(records, nhdu))              # update FITS_header object    
    push!(FITS_data, _cast_data(nhdu, hdutype, data))             # update FITS_data object 
    
    FITS = [FITS_HDU(filename, i, FITS_headers[i], FITS_data[i]) for i=1:nhdu]
    
    _fits_save(FITS)
       
    return FITS
    
end


"""
    fits_copy(filenameA [, filenameB="" [; protect=true]])

Copy "filenameA" to "filenameB" (with mandatory ".fits" extension)
#### Examples:
```

fits_copy("T01.fits")
'T01.fits' was saved as 'T01 - Copy.fits'

fits_copy("T01.fits", "T01a.fits")
filename (T01a.fits) in use (overwrite protected)

fits_copy("T01.fits", "T01a.fits"; protect=false)
'T01.fits' was saved as 'T01a.fits'

```
"""
function fits_copy(filenameA::String, filenameB::String=" "; protect=true)
    
    o = _fits_read_IO(filenameA)   
    f = cast_FITS_name(filenameA)
    
    filenameB = filenameB == " " ? "$(f.name) - Copy.fits" : filenameB
    
    _validate_FITS_name(filenameB)
    
    strA = "'$filenameA' was saved as '$filenameB'"
    strB = "'$filenameA': copy failed"
    
    return _isavailable(filenameB, protect) ? (_fits_write_IO(o,filenameB); strA) : strB
    
end


"""
    fits_add_key(filename, hduindex, key, value, comment)

Add a header record of given 'key, value and comment' to 'HDU[hduindex]' of file with name 'filename'
#### Example:
```

filnam="minimal.fits"
fits_create(filnam;protect=false)
fits_add_key(filnam, 1, "EXTEND2", true, "FITS dataset may contain extension")
 'EXTEND2': key added; new record: 'EXTEND2 =                    T / FITS dataset may contain extension             '

f = fits_read(filnam)
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
 EXTEND2 =                    T / FITS dataset may contain extension             
 END

```
"""
function fits_add_key(filename::String, hduindex::Int, key::String, val::Real, com::String)
    
    o = _fits_read_IO(filename)  
    
    nhdu = _hdu_count(o)
    
    FITS_headers = [_read_header(o,i) for i=1:nhdu] 
       FITS_data = [_read_data(o,i) for i=1:nhdu] 
    
    str = "'$key' key truncated at 8 characters (FITS standard)"
    length(key) < 9 ? key = Base.strip(key) : (println(str); key = key[1:8])
    val = val == true ? "T" : val == false ? "F" : val
     
    H = FITS_headers[hduindex]
        Base.haskey(H.maps,key) && return println("'$key': key in use (use different key name or edit key)")
    i = Base.get(H.maps, "END", "Error: key not found")
        H.records[i] = rpad(key,8) * "= " * lpad(val,20) * " / " * rpad(com,47)
        Base.push!(H.records, "END" * Base.repeat(" ",77))
    H = _cast_header(H.records, hduindex)                               # here we cast FITS_headers[hduindex] 
    
    _fits_save([FITS_HDU(filename, i, FITS_headers[i], FITS_data[i]) for i ∈ eachindex(FITS_headers)])
    
    return println("'$key': key added; new record: '$(H.records[i])'")
    
end


"""
    fits_edit_key(filename, hduindex, key, value, comment)

Edit a header record of given 'key, value and comment' to 'HDU[hduindex]' of file with name 'filename'
#### Example:
```

filnam="minimal.fits"
fits_create(filnam;protect=false)
fits_add_key(filnam, 1, "EXTEND12", true, "FITS dataset may contain extension")
'EXTEND12': key added; new record: 'EXTEND12=                    T / FITS dataset may contain extension             '

fits_edit_key(filnam, 1, "EXTEND12", true, "FITS dataset new comment")
'EXTEND12': key edited; new record: 'EXTEND12=                    F / FITS dataset new comment                       '

```
"""
function fits_edit_key(filename::String, hduindex::Int, key::String, val::Real, com::String)
    
    o = _fits_read_IO(filename)   
    
    nhdu = _hdu_count(o)
    
    FITS_headers = [_read_header(o,i) for i=1:nhdu] 
       FITS_data = [_read_data(o,i) for i=1:nhdu]  
    
    key = strip(key)
    val = val == true ? "T" : val == false ? "F" : val
    res = ["SIMPLE","BITPIX","NAXIS","NAXIS1","NAXIS2","NAXIS3","BZERO","END"]
    
    key ∈ res && return println("'$key': cannot be edited (key protected under FITS standard)")
     
    H = FITS_headers[hduindex]
        Base.haskey(H.maps, key) ||  return println("'$key': cannot be deleted (key not found)")
    i = Base.get(H.maps, key, "error: key not found")
        H.records[i] = Base.rpad(key,8) * "= " * Base.lpad(val,20) * " / " * Base.rpad(com,47)
    H = _cast_header(H.records, hduindex) 
    
    _fits_save([FITS_HDU(filename, i, FITS_headers[i], FITS_data[i]) for i ∈ eachindex(FITS_headers)])
    
    return println("'$key': key edited; new record: '$(H.records[i])'")
    
end


"""
    fits_delete_key(filename, hduindex, key)

Delete a header record of given `key`, `value` and `comment` to `FITS_HDU[hduindex]` of file with name  'filename'
#### Examples:

```
filnam="minimal.fits"
fits_create(filnam;protect=false)
fits_add_key(filnam, 1, "EXTEND123", true, "FITS dataset may contain extension")
 'EXTEND123' key truncated at 8 characters (FITS standard)
 'EXTEND12': key added; new record: 'EXTEND12=                    T / FITS dataset may contain extension             '

fits_add_key(filnam, 1, "EXTEND12", true, "FITS dataset may contain extension")
 'EXTEND12': key in use (use different key name or edit key)

fits_delete_key(filnam, 1, "EXTEND12")
 'EXTEND12': key deleted

fits_delete_key(filnam, 1, "EXTEND12")
 'EXTEND12': cannot be deleted (key not found)

fits_delete_key(filnam, 1, "NAXIS")
 'NAXIS': cannot be deleted (key protected under FITS standard)

```
"""
function fits_delete_key(filename::String, hduindex::Int, key::String)
    
    o = _fits_read_IO(filename)  
    
    nhdu = _hdu_count(o)
    
    FITS_headers = [_read_header(o,i) for i=1:nhdu] 
       FITS_data = [_read_data(o,i) for i=1:nhdu]   
    
    key = strip(key)
    res = ["SIMPLE","BITPIX","NAXIS","NAXIS1","NAXIS2","NAXIS3","BZERO","BSCALE","END"]
    
    key ∈ res && return println("'$key': cannot be deleted (key protected under FITS standard)")
    
    H = FITS_headers[hduindex]
        Base.haskey(H.maps, key) || return println("'$key': cannot be deleted (key not found)") 
    i = Base.get(H.maps, key, "error: key not found")
        Base.splice!(H.records,i)
    H = _cast_header(H.records, hduindex)                           # here we cast FITS_headers[hduindex] 
    
    _fits_save([FITS_HDU(filename, i, FITS_headers[i], FITS_data[i]) for i ∈ eachindex(FITS_headers)])
    
    return println("'$key': key deleted")
    
end
