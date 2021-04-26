# ......................................... FITS public sector .................................................................
# .................................................... fits_create ...................................................
"""
    fits_create(filename [, data [; protect=true]])

Create FITS file of given filename [, optional data block [, default overwrite protection]] and return Array of HDUs.
Key:
* `protect::Bool`: overwrite protection
#### Examples:
```
strExample = "minimal.fits"
fits_create(strExample;protect=false)

f = fits_read(strExample)
a = f[1].dataobject.data
b = f[1].header.keys
println(a);println(b)
  Any[]
  ["SIMPLE", "NAXIS", "EXTEND", "COMMENT", "END"]

strExample = "example.fits"
data = [0x0000043e, 0x0000040c, 0x0000041f]
fits_create(strExample, data; protect=false)

f = fits_read(strExample)
a = f[1].dataobject.data
b = f[1].header.keys
println(a);println(b)
  UInt32[0x0000043e, 0x0000040c, 0x0000041f]
  ["SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "BZERO", "BSCALE", "EXTEND", "COMMENT", "END"]
```
"""
function fits_create(filename::String, data=[]; protect=true)

    strErr = "FitsError: '$filename': creation failed (filename in use - set ';protect=false' to overrule overwrite protection)"

    _validate_FITS_name(filename)
    _isavailable(filename, protect) || error(strErr)

    nhdu = 1
    hdutype = "PRIMARY"

       FITS_data = [_cast_data(i, hdutype, data) for i=1:nhdu]
    FITS_headers = [_cast_header(_PRIMARY_input(FITS_data[i]), i) for i=1:nhdu]

    FITS = [FITS_HDU(filename, i, FITS_headers[i], FITS_data[i]) for i=1:nhdu]

    return _fits_save(FITS)

end
# test ...
function fits_create()

    strExample = "minimal.fits"
    fits_create(strExample, protect=false)

    f = fits_read(strExample)
    a = f[1].header.keys[1]  == "SIMPLE"
    b = f[1].dataobject.data == Any[]
    c = get(Dict(f[1].header.dict),"SIMPLE",0)
    d = get(Dict(f[1].header.dict),"NAXIS",0) == 0;

    rm(strExample)

    test = .![a, b, c, d];

    return !convert(Bool,sum(test))

end


# .................................................... fits_info ...................................................
"""
    fits_info(hdu)

Print metafinformation of given `FITS_HDU`
#### Example:
```
strExample = "minimal.fits"
fits_create(strExample)

f = fits_read(strExample)
p = fits_info(f[1]); println(p)

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
function fits_info(hdu::FITS_HDU)

    typeof(hdu) <: FITS_HDU || error("FitsWarning: FITS_HDU not found")

    info = [
        "\r\nFile: " * hdu.filename,
        "hdu: " * Base.string(hdu.hduindex),
        "hdutype: " * hdu.dataobject.hdutype,
        "DataType: " * Base.string(Base.eltype(hdu.dataobject.data)),
        "Datasize: " * Base.string(Base.size(hdu.dataobject.data)),
        "\r\nMetainformation:"
        ]

    records = hdu.header.records

    Base.append!(info, records)

    return Base.join(info .* "\r\n")

end
#  test ...
function fits_info()

    strExample = "minimal.fits"
    fits_create(strExample; protect=false)

    f = fits_read(strExample)
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

    test = fits_info(f[1]) == Base.join(info .* "\r\n")

    rm(strExample)

    return test

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
    fits_create(strExample, protect=false)

    f = fits_read(strExample)
    a = f[1].header.keys[1]  == "SIMPLE"
    b = f[1].dataobject.data == Any[]
    c = get(Dict(f[1].header.dict),"SIMPLE",0)
    d = get(Dict(f[1].header.dict),"NAXIS",0) == 0;

    rm(strExample)

    test = .![a, b, c, d];

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
fits_create(strExample, data, protect=false)

f = fits_read(strExample)
a = Float16[1.01E-6,2.0E-6,3.0E-6,4.0E-6,5.0E-6]
b = [0x0000043e, 0x0000040c, 0x0000041f, 0x0000042e, 0x0000042f]
c = [1.23,2.12,3.,4.,5.]
d = ['a','b','c','d','e']
e = ["a","bb","ccc","dddd","ABCeeaeeEEEEEEEEEEEE"]
data = [a,b,c,d,e]
fits_extend(strExample, data, "TABLE")

f = fits_read(strExample)
f[2].dataobject.data
  5-element Vector{String}:
   "1.0e-6 1086 1.23 a a                    "
   "2.0e-6 1036 2.12 b bb                   "
   "3.0e-6 1055 3.0  c ccc                  "
   "4.0e-6 1070 4.0  d dddd                 "
   "5.0e-6 1071 5.0  e ABCeeaeeEEEEEEEEEEEE "

rm(strExample); f = data = a = b = c = d = e = nothing
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
    fits_create(strExample, data, protect=false)

    f = fits_read(strExample)
    a = Float16[1.01E-6,2.0E-6,3.0E-6,4.0E-6,5.0E-6]
    b = [0x0000043e, 0x0000040c, 0x0000041f, 0x0000042e, 0x0000042f]
    c = [1.23,2.12,3.,4.,5.]
    d = ['a','b','c','d','e']
    e = ["a","bb","ccc","dddd","ABCeeaeeEEEEEEEEEEEE"]
    data = [a,b,c,d,e]
    fits_extend(strExample, data, "TABLE")

    f = fits_read(strExample)
    a = f[1].header.keys[1]  == "SIMPLE"
    b = f[1].dataobject.data[1] == 0x0000043e
    c = f[2].header.keys[1]  == "XTENSION"
    d = f[2].dataobject.data[1] == "1.0e-6 1086 1.23 a a                    "
    e = get(Dict(f[2].header.dict),"NAXIS",0) == 2

    rm(strExample)

    test = .![a, b, c, d, e]

    return !convert(Bool,sum(test))

end

# .................................................... fits_copy ...................................................

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

# .................................................... fits_add_key ...................................................

"""
    fits_add_key(filename, hduindex, key, value, comment)

Add a header record of given 'key, value and comment' to 'HDU[hduindex]' of file with name 'filename'
#### Example:
```
strExample="minimal.fits"
fits_create(strExample;protect=false)
fits_add_key(strExample, 1, "KEYNEW1", true, "FITS dataset may contain extension")

f = fits_read(strExample)
p = fits_info(f[1]); println(p)

  File: minimal.fits
  HDU: 1
  DataType: Any
  Datasize: (0,)

  Metainformation:
  SIMPLE  =                    T / file does conform to FITS standard
  NAXIS   =                    0 / number of data axes
  EXTEND  =                    T / FITS dataset may contain extensions
  COMMENT    Basic FITS file     / http://fits.gsfc.nasa.gov/iaufwg
  KEYNEW1 =                    T / FITS dataset may contain extension
  END
```
"""
function fits_add_key(filename::String, hduindex::Int, key::String, val::Any, com::String)

    o = _fits_read_IO(filename)

    nhdu = _hdu_count(o)

    FITS_headers = [_read_header(o,i) for i=1:nhdu]
       FITS_data = [_read_data(o,i) for i=1:nhdu]

    key = _format_key(key)

    h = FITS_headers[hduindex]
    Base.get(h.maps, key, 0) > 0 && return println("FitsError: '$key': key in use (use different name or edit key)")

    newrecords = _fits_new_records(key, val, com)

    Base.pop!(h.records)
   [Base.push!(h.records, newrecords[i]) for i ∈ eachindex(newrecords)]
    Base.push!(h.records, "END" * Base.repeat(" ",77))

    FITS_headers[hduindex] = _cast_header(h.records, hduindex)

    FITS = [FITS_HDU(filename, i, FITS_headers[i], FITS_data[i]) for i=1:nhdu]

    return _fits_save(FITS)

end
# test ...
function fits_add_key()

    strExample="minimal.fits"
    fits_create(strExample;protect=false)
    fits_add_key(strExample, 1, "KEYNEW1", true, "FITS dataset may contain extension")

    f = fits_read(strExample)
    i = get(f[1].header.maps,"KEYNEW1",0)
    r = f[1].header.records;

    test = r[i] == "KEYNEW1 =                    T / FITS dataset may contain extension             "

    rm(strExample)

    return test

end

"""
    fits_edit_key(filename, hduindex, key, value, comment)

Edit a header record of given 'key, value and comment' to 'HDU[hduindex]' of file with name 'filename'
#### Example:
```
    strExample="minimal.fits"
    fits_create(strExample;protect=false)
    fits_add_key(strExample, 1, "KEYNEW1", true, "FITS dataset may contain extension")
    fits_edit_key(strExample, 1, "KEYNEW1", false, "comment has changed")

    f = fits_read(strExample)
    i = get(f[1].header.maps,"KEYNEW1",0)
    f[1].header.records[i]
      "KEYNEW1 =                    F / comment has changed                            "

    rm(strExample)
```
"""
function fits_edit_key(filename::String, hduindex::Int, key::String, val::Any, com::String)

    o = _fits_read_IO(filename)

    nhdu = _hdu_count(o)

    FITS_headers = [_read_header(o,i) for i=1:nhdu]
       FITS_data = [_read_data(o,i) for i=1:nhdu]

    key = _format_key(key)
    res = ["SIMPLE","BITPIX","NAXIS","NAXIS1","NAXIS2","NAXIS3","BZERO","END"]
    key ∈ res && return println("FitsError: '$key': cannot be edited (key protected under FITS standard)")

    h = FITS_headers[hduindex]
    i = Base.get(h.maps, key, 0)
    i == 0 && return println("FitsError: '$key': key not found")

    nold = length(h.keys)
    nobs = length(_fits_obsolete_records(h,i))
    newrecords = _fits_new_records(key, val, com)
    oldrecords = h.records[i+nobs:end]
   [Base.pop!(h.records) for j=i:nold]
   [Base.push!(h.records, newrecords[i]) for i ∈ eachindex(newrecords)]
   [Base.push!(h.records, oldrecords[i]) for i ∈ eachindex(oldrecords)]

    FITS_headers[hduindex] = _cast_header(h.records, hduindex)

    FITS = [FITS_HDU(filename, i, FITS_headers[i], FITS_data[i]) for i=1:nhdu]

    return _fits_save(FITS)

end
# test ...
function fits_edit_key()

    strExample="minimal.fits"
    fits_create(strExample;protect=false)
    fits_add_key(strExample, 1, "KEYNEW1", true, "FITS dataset may contain extension")
    fits_edit_key(strExample, 1, "KEYNEW1", false, "comment has changed")

    f = fits_read(strExample)
    i = get(f[1].header.maps,"KEYNEW1",0)
    r = f[1].header.records;

    test = r[i] == "KEYNEW1 =                    F / comment has changed                            "

    rm(strExample)

    return test

end

"""
    fits_delete_key(filename, hduindex, key)

Delete a header record of given `key`, `value` and `comment` to `FITS_HDU[hduindex]` of file with name  'filename'
#### Examples:
```
strExample="minimal.fits"
fits_create(strExample;protect=false)
fits_add_key(strExample, 1, "KEYNEW1", true, "this is record 5")

f = fits_read(strExample)
get(f[1].header.maps,"KEYNEW1",0)
  5

fits_delete_key(strExample, 1, "KEYNEW1")

f = fits_read(strExample)
get(f[1].header.maps,"KEYNEW1",0)
  0

fits_delete_key(filnam, 1, "NAXIS")
 'NAXIS': cannot be deleted (key protected under FITS standard)
```
"""
function fits_delete_key(filename::String, hduindex::Int, key::String)

    o = _fits_read_IO(filename)

    nhdu = _hdu_count(o)

    FITS_headers = [_read_header(o,i) for i=1:nhdu]
       FITS_data = [_read_data(o,i) for i=1:nhdu]

    key = _format_key(key)
    res = ["SIMPLE","BITPIX","NAXIS","NAXIS1","NAXIS2","NAXIS3","BZERO","END"]
    key ∈ res && return println("FitsError: '$key': cannot be edited (key protected under FITS standard)")

    h = FITS_headers[hduindex]
    i = Base.get(h.maps, key, 0)
    i == 0 && return println("FitsError: '$key': key not found")

    nold = length(h.keys)
    nobs = length(_fits_obsolete_records(h,i))
    oldrecords = h.records[i+nobs:end]
   [Base.pop!(h.records) for j=i:nold]
   [Base.push!(h.records, oldrecords[i]) for i ∈ eachindex(oldrecords)]

    FITS_headers[hduindex] = _cast_header(h.records, hduindex)

    FITS = [FITS_HDU(filename, i, FITS_headers[i], FITS_data[i]) for i=1:nhdu]

    _fits_save(FITS)

    return FITS

end
# test ...
function fits_delete_key()

    strExample="minimal.fits"
    fits_create(strExample;protect=false)
    fits_add_key(strExample, 1, "KEYNEW1", true, "FITS dataset may contain extension")

    f = fits_read(strExample)
    i = get(f[1].header.maps,"KEYNEW1",0)

    test1 = i == 5

    fits_delete_key(strExample, 1, "KEYNEW1")

    f = fits_read(strExample)
    i = get(f[1].header.maps,"KEYNEW1",0)

    test2 = i == 0

    test = .![test1, test2]

    rm(strExample)

    return !convert(Bool,sum(test))

end

"""
    fits_rename_key(filename, hduindex, keyold, kewnew)

Rename the key of a header record of file with name 'filename'
#### Example:
```
strExample="minimal.fits"
fits_create(strExample;protect=false)
fits_add_key(strExample, 1, "KEYNEW1", true, "this is record 5")

f = fits_read(strExample)
get(f[1].header.maps,"KEYNEW1",0)
  5

fits_rename_key(strExample, 1,"KEYNEW1", "KEYNEW2")

f = fits_read(strExample)
get(f[1].header.maps,"KEYNEW2",0)
  5

rm(strExample)
```
"""
function fits_rename_key(filename::String, hduindex::Int, keyold::String, keynew::String)

    o = _fits_read_IO(filename)

    nhdu = _hdu_count(o)

    FITS_headers = [_read_header(o,i) for i=1:nhdu]
       FITS_data = [_read_data(o,i) for i=1:nhdu]

    keyold = _format_key(keyold)
       res = ["SIMPLE","BITPIX","NAXIS","NAXIS1","NAXIS2","NAXIS3","BZERO","END"]
    keyold ∈ res && return println("FitsWarning: '$keyold': cannot be renamed (key protected under FITS standard)")

    h = FITS_headers[hduindex]
    i = Base.get(h.maps, keyold, 0)
    i == 0 && return println("FitsError: '$keyold': key not found")
    Base.get(h.maps, keynew, 0) > 0 && return println("FitsWarning: '$keynew': key in use (use different name or edit key)")

    h.records[i] = rpad(keynew,8) * h.records[i][9:80]

    FITS_headers[hduindex] = _cast_header(h.records, hduindex)

    FITS = [FITS_HDU(filename, i, FITS_headers[i], FITS_data[i]) for i=1:nhdu]

    _fits_save(FITS)

    return FITS

end
# test ...
function fits_rename_key()

    strExample="minimal.fits"
    fits_create(strExample;protect=false)
    fits_add_key(strExample, 1, "KEYNEW1", true, "this is record 5")

    f = fits_read(strExample)
    i = get(f[1].header.maps,"KEYNEW1",0)

    test1 = i == 5

    fits_rename_key(strExample, 1,"KEYNEW1", "KEYNEW2")

    f = fits_read(strExample)
    i = get(f[1].header.maps,"KEYNEW2",0)

    test2 = i == 5

    test = .![test1, test2]

    rm(strExample)

    return !convert(Bool,sum(test))

end
