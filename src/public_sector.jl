# ......................................... FITS public sector .................................................................

# .................................................... fits_info ...................................................
"""
    fits_info(FITS_HDU; printformat=true)

Print metafinformation of given `FITS_HDU`

Key:
* `printformat::Bool`: output formatted by function `print`
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
function fits_info(FITS_HDU; printformat=true)

    info = [
        "\r\nFile: " * FITS_HDU.filename,
        "hdu: " * Base.string(FITS_HDU.hduindex),
        "hdutype: " * FITS_HDU.dataobject.hdutype,
        "DataType: " * Base.string(Base.eltype(FITS_HDU.dataobject.data)),
        "Datasize: " * Base.string(Base.size(FITS_HDU.dataobject.data)),
        "\r\nMetainformation:"
        ]

    records = FITS_HDU.header.records

    Base.append!(info,_rm_blanks(records))

    return printformat ? print(Base.join(info .* "\r\n")) : Base.join(info .* "\r\n")

end
#  test ...
function fits_info()

    strExample = "minimal.fits"
    f = fits_create(strExample, protect=false)

    info = [
            "\r\nFile: minimal.fits",
            "hdu: 1",
            "hdutype: PRIMARY",
            "DataType: Any",
            "Datasize: (0,)",
            "\r\nMetainformation:",
            "SIMPLE  =                    T / file does conform to FITS standard             ",
            "NAXIS   =                    0 / number of data axes                            ",
            "EXTEND  =                    T / FITS dataset may contain extensions            ",
            "COMMENT    Primary FITS HDU    / http://fits.gsfc.nasa.gov/iaufwg               ",
            "END                                                                             "
            ]

    test = fits_info(f[1]; printformat=false) == Base.join(info .* "\r\n")

    rm(strExample)

    return test

end

# .................................................... fits_create ...................................................
"""
    fits_create(filename [, data [; protect=true]])

Create FITS file of given filename [, optional data block [, default overwrite protection]] and return Array of HDUs.
Key:
* `protect::Bool`: overwrite protection
#### Examples:
```
strExample = "minimal.fits"
f = fits_create(strExample;protect=false)
a = f[1].dataobject.data
b = f[1].header.keys
println(a);println(b)
  Any[]
  ["SIMPLE", "NAXIS", "EXTEND", "COMMENT", "END"]

strExample = "example.fits"
data = [0x0000043e, 0x0000040c, 0x0000041f]
f = fits_create(strExample, data; protect=false)
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
# test ...
function fits_create()

    strExample = "minimal.fits"
    f = fits_create(strExample, protect=false)

    test1 = f[1].header.keys[1]  == "SIMPLE"
    test2 = f[1].dataobject.data == Any[]
    test3 = get(Dict(f[1].header.dict),"SIMPLE",0)
    test4 = get(Dict(f[1].header.dict),"NAXIS",0) == 0

    rm(strExample)

    test = .![test1, test2, test3, test4]

    return !convert(Bool,sum(test))

end

# .................................................... fits_read ...................................................
"""
    fits_read(filename)

Read FITS file and return Array of `FITS_HDU`s
#### Example:
```
strExample = "minimal.fits"
fits_create(strExample;protect=false)
f = fits_read(strExample)
f[1].dataobject.data
  Any[]

rm(strExample); f = nothing
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
# test ...
function fits_read()

    strExample = "minimal.fits"
    f = fits_create(strExample, protect=false)
    f = nothing
    f = fits_read(strExample)

    test1 = f[1].header.keys[1]  == "SIMPLE"
    test2 = f[1].dataobject.data == Any[]
    test3 = get(Dict(f[1].header.dict),"SIMPLE",0)
    test4 = get(Dict(f[1].header.dict),"NAXIS",0) == 0

    rm(strExample); f = nothing

    test = .![test1, test2, test3, test4]

    return !convert(Bool,sum(test))

end


# .................................................... fits_extend ...................................................
"""
    fits_extend(filename, data_extend, hdutype="IMAGE")

Extend the FITS file of given filename with the data of `hdutype` from `data_extend`  and return Array of HDUs.
#### Examples:
```
strExample = "test_example.fits"
data = [0x0000043e, 0x0000040c, 0x0000041f]
f = fits_create(strExample, data, protect=false)
table1 = Float16[1.01E-6,2.0E-6,3.0E-6,4.0E-6,5.0E-6]
table2 = [0x0000043e, 0x0000040c, 0x0000041f, 0x0000042e, 0x0000042f]
table3 = [1.23,2.12,3.,4.,5.]
table4 = ['a','b','c','d','e']
table5 = ["a","bb","ccc","dddd","ABCeeaeeEEEEEEEEEEEE"]
data = [table1,table2,table3,table4,table5]
f = fits_extend(strExample, data, "TABLE")
f[2].dataobject.data
  5-element Vector{String}:
   "1.0e-6 1086 1.23 a a                    "
   "2.0e-6 1036 2.12 b bb                   "
   "3.0e-6 1055 3.0  c ccc                  "
   "4.0e-6 1070 4.0  d dddd                 "
   "5.0e-6 1071 5.0  e ABCeeaeeEEEEEEEEEEEE "

rm(strExample); f = nothing; data =nothing
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
# test ...
function fits_extend()

    strExample = "test_example.fits"
    data = [0x0000043e, 0x0000040c, 0x0000041f]
    f = fits_create(strExample, data, protect=false)

    table1 = Float16[1.01E-6,2.0E-6,3.0E-6,4.0E-6,5.0E-6]
    table2 = [0x0000043e, 0x0000040c, 0x0000041f, 0x0000042e, 0x0000042f]
    table3 = [1.23,2.12,3.,4.,5.]
    table4 = ['a','b','c','d','e']
    table5 = ["a","bb","ccc","dddd","ABCeeaeeEEEEEEEEEEEE"]
    data = [table1,table2,table3,table4,table5]
    f = fits_extend(strExample, data, "TABLE")

    test1 = f[1].header.keys[1]  == "SIMPLE"
    test2 = f[1].dataobject.data[1] == 0x0000043e
    test3 = f[2].header.keys[1]  == "XTENSION"
    test4 = f[2].dataobject.data[1] == "1.0e-6 1086 1.23 a a                    "
    test5 = get(Dict(f[2].header.dict),"NAXIS",0) == 2

    rm(strExample); f = nothing; data = nothing

    test = .![test1, test2, test3, test4, test5]

    return !convert(Bool,sum(test))

end


"""
    fits_copy(filenameA [, filenameB="" [; protect=true]])

Copy "filenameA" to "filenameB" (with mandatory ".fits" extension)
Key:
* `protect::Bool`: overwrite protection
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

    _validate_FITSname(filenameB)

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
function fits_add_key(filename::String, hduindex::Int, key::String, val::Union{Real,String,Char}, com::String)

    o = _fits_read_IO(filename)

    nhdu = _hdu_count(o)

    FITS_headers = [_read_header(o,i) for i=1:nhdu]
       FITS_data = [_read_data(o,i) for i=1:nhdu]

    newrecords = _fits_new_records(key, val, com)
    nrec = length(newrecords)

    H = FITS_headers[hduindex]
        Base.haskey(H.maps,key) && return println("'$key': key in use (use different key name or edit key)")
    i = Base.get(H.maps, "END", 0)
        i > 0 ? H.records[i] = newrecords[1] : error("Error: key not found")
        nrec > 1 ? [Base.push!(H.records, newrecords[i]) for i=2:nrec] : 0
        Base.push!(H.records, "END" * Base.repeat(" ",77))

    FITS_headers[hduindex] = _cast_header(H.records, hduindex)

    FITS = [FITS_HDU(filename, i, FITS_headers[i], FITS_data[i]) for i=1:nhdu]

    _fits_save(FITS)

    return FITS

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
function fits_edit_key(filename::String, hduindex::Int, key::String, val::Any, com::String)
println("hoi")
    o = _fits_read_IO(filename)

    nhdu = _hdu_count(o)

    FITS_headers = [_read_header(o,i) for i=1:nhdu]
       FITS_data = [_read_data(o,i) for i=1:nhdu]

    key = string(strip(key))
    res = ["SIMPLE","BITPIX","NAXIS","NAXIS1","NAXIS2","NAXIS3","BZERO","END"]
    key ∈ res && return println("'$key': cannot be edited (key protected under FITS standard)")

    val = _is_recordvalue_charstring(val::Any) ? val : typeof(val) <: Real ? val : error("FitsError: "not a valid value type")

    newrecords = _fits_new_records(key, val, com)
    nrec = length(newrecords)

    H = FITS_headers[hduindex]
        Base.haskey(H.maps, key) ||  error("FitsError: '$key': cannot be edited (key not found)")
    i = Base.get(H.maps, key, 0)
        H.records[i] = newrecords[1]
        nrec > 1 ? [Base.push!(H.records, newrecords[i]) for i=2:nrec] : 0

    FITS_headers[hduindex] = _cast_header(H.records, hduindex)

    FITS = [FITS_HDU(filename, i, FITS_headers[i], FITS_data[i]) for i=1:nhdu]

    _fits_save(FITS)

    return FITS

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

    key = join(strip(key))
    res = ["SIMPLE","BITPIX","NAXIS","NAXIS1","NAXIS2","NAXIS3","BZERO","BSCALE","END"]

    key ∈ res && return println("'$key': cannot be deleted (key protected under FITS standard)")

    H = FITS_headers[hduindex]
        Base.haskey(H.maps, key) || error("FitsError: '$key': cannot be deleted (key not found)")
    i = Base.get(H.maps, key, 0)
        Base.splice!(H.records,i)

    FITS_headers[hduindex] = _cast_header(H.records, hduindex)

    FITS = [FITS_HDU(filename, i, FITS_headers[i], FITS_data[i]) for i=1:nhdu]

    _fits_save(FITS)

    return FITS


end
