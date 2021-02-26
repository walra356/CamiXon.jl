using FITSIO

# A FITS file consists of a name and one or more header-data-units (HDUs) concatenated one after the other. The FITS object is a collection of these HDUs.

# The file 'filnamA.fits' is opened for 'reading only' by the command fileA = FITS(filnamA,'r').
# The second argument can be 'r' (read-only; default), 'r+' (read-write) or 'w' (read/create-write). 
# In the 'write' mode, a  new file is created or any existing file of the same name is overwritten. 

# The file fileA consists of one or more HDU's, fileA[1](, fileA[2],...). 
# HDUs consist of data blocks with a header containing metainformation, accessible by the commands read() and read_header().
# In the given example this becomes: dataA = read(fileA[1]) and metaInfoA = read_header(fileA[1])"
# The 3D array dataA represents a stack of images. The first image is obtained by the command imgA = dataA[:, :,1]. 

function _file_exists(filnam::String)
    
    return Base.Filesystem.isfile(filnam) ? true : false
    
end

function _key_size_ok(key::String)
    
    return length(key) < 9 ? true : false
    
end 

function _key_comment_ok(comment::String)
    
    return length(comment) < 48 ? true : false
    
end       
    
function _key_exists(filnam::String,key::String)
    
    f = FITSIO.FITS(filnam,"r+")       # open file (read-write mode)
    m = FITSIO.read_header(f[1])       # read hdu header from file
    
    Base.close(f)
        
    return Base.haskey(m, key) ? true : false
    
end         
    
function _key_available(filnam::String,key::String)
    
    f = FITSIO.FITS(filnam,"r+")       # open file (read-write mode)
    m = FITSIO.read_header(f[1])       # read hdu header from file
    
    Base.close(f)
        
    return !Base.haskey(m, key) ? true : false
    
end   
    
function _key_reserved(filnam::String,key::String)
    
    f = FITSIO.FITS(filnam,"r+")       # open file (read-write mode)
    m = FITSIO.read_header(f[1])       # read hdu header from file
    r = FITSIO.reserved_key_indices(m) 
    i = m.map[key]
    
    Base.close(f)
    
    return i ∈ r ? true : false
    
end

"""
    decompose_filnam(str)

Decompose filename into its name (and, if present, extension, prefix and numerator).
#### Examples:
```
strExample = "T23.01.fits"

dict = decompose_filnam(strExample)
Dict{String,String} with 4 entries:
  "Extension" => ".fits"
  "Numerator" => "01"
  "Prefix"    => "T23."
  "Name"      => "T23.01"

get(dict,"Numerator","Numerator: Key absent")
"01"

get(dict,"Wild","Key absent")
"Key absent"
```
"""
function decompose_filnam(str::String)
    
    ne = Base.findlast('.',str)                # ne: first digit of extension
    nl = length(str)                           # ne: length of file name including extension
    ne == nothing ? extension = false : extension = true

    if extension
        strNam = str[1:ne-1]
        strExt = str[ne:nl]
        Base.Unicode.isdigit(str[ne-1]) ? n = ne-1 : n = nothing
        o = [("Name", strNam), ("Extension", strExt)]
    else
        strNam = str[1:ne-1]
        strExt = nothing
        Base.Unicode.isdigit(str[nl]) ? n = nl : n = nothing
        o = [("Name", strNam)]
    end

    if n != nothing
        strNum = ""
        while Base.Unicode.isdigit(str[n])
            strNum = str[n] * strNum
            n -= 1
        end
        strPre = str[1:n]
        Base.append!(o,[("Prefix", strPre),("Numerator", strNum)])
    end

    return Dict(o)

end

"""
    fits_combine(str1, str2 [; info=false])

Combine a series of .fits files into a single .fits file.
#### Example:
```
fits_combine("T01.fits", "T22.fits"; info=false)
T01-T22.FITS: file was created (for more information set info=true)
```
"""
function fits_combine(filnamFirst::String, filnamLast::String; info=false)
    
    _file_exists(filnamFirst) || return "$filnamFirst: file not found"
    _file_exists(filnamLast) || return "$filnamLast: file not found"

    filnamFirst = uppercase(filnamFirst)
    filnamLast = uppercase(filnamLast)

    d = decompose_filnam(filnamFirst)
    strPre = get(d,"Prefix","Error: no prefix")
    strNum = get(d,"Numerator","Error: no Numerator")
    strExt = get(d,"Extension","Error: no extension")
    valNum = parse(Int,strNum )
    numLeadingZeros = length(strNum) - length(string(valNum))

    d = decompose_filnam(filnamLast)
    strPre2 = get(d,"Prefix","Error: no prefix")
    strNum2 = get(d,"Numerator","Error: no Numerator")
    strExt2 = get(d,"Extension","Error: no extension")
    valNum2 = parse(Int,strNum2 )
    numLeadingZeros2 = length(strNum2) - length(string(valNum2))

    if strPre ≠ strPre2
        error(strPre * " ≠ " * strPre2 * " (prefixes must be identical)")
    elseif strExt ≠ strExt2
        error(strExt * " ≠ " * strExt2 * " (file extensions must be identical)")
    elseif strExt ≠ ".FITS"
        error("file extension must be '.fits'")
    end

    numFiles = 1 + valNum2 - valNum
    fileFirst = FITSIO.FITS(filnamFirst)
    metaInfo = FITSIO.read_header(fileFirst[1])
    dataFirst = Base.read(fileFirst[1])  # read an image from disk
    Base.close(fileFirst)
    t = typeof(dataFirst[1,1,1])
    s = size(dataFirst)
    dataStack =  Array{t,3}(undef, s[1], s[2] , numFiles)

    itr = valNum:valNum2
    filnamNext = filnamFirst
    for i ∈ itr
        l = length(filnamNext)
        filnamNext = strPre * "0"^numLeadingZeros * string(i) * ".fits"
        if l < length(filnamNext)
            numLeadingZeros = numLeadingZeros -1
            filnamNext = strPre * "0"^numLeadingZeros * string(i) * ".fits"
        end
        fileNext = FITS(filnamNext)
        dataNext = Base.read(fileNext[1])  # read an image from disk
        Base.close(fileNext)
        dataStack[:, :,i] = dataNext[:, :,1]
    end

    filnamOut = strPre * strNum * "-" * strPre * strNum2 * strExt
    fileOut = FITSIO.FITS(filnamOut,"w")
    Base.write(fileOut, dataStack; header=metaInfo)
    if info
        println("Output fileOut:\r\n", fileOut)
        println("\r\nOutput fileOut[1]:\r\n", fileOut[1])
        Base.close(fileOut)
        println("\r\nmetaInformation:\r\n", metaInfo)
    else
        Base.close(fileOut)
        return filnamOut * ": file was created (for more information set info=true)"
    end
end


"""
    fits_copy(filnam [, filnam2="" [; protect=true]])

Copy "filnam.fits" to "filnam2.fits"
#### Examples:
```
fits_copy("T01.fits")
T01.fits was saved as T01 - Copy.fits

fits_copy("T01.fits", "T01a.fits")
"T01a.fits: filname in use (overwrite protected)"

fits_copy("T01.fits", "T01a.fits"; protect=false)
T01.fits was saved as T01a.fits
```
"""
function fits_copy(filnam, filnam2=""; protect=true)
    
    _file_exists(filnam) || return "$filnam: file not found"
    
    _file_exists(filnam2) &&  protect && return "$filnam2: filname in use (overwrite protected)"

    d1 = decompose_filnam(filnam)
    
    strNam = Base.get(d1,"Name","Error: no name")
    strExt = Base.get(d1,"Extension","Error: no extension")
        
    if filnam2 == ""
        filnam2 = strNam * " - Copy" * strExt
    else
        d2 = decompose_filnam(filnam2)
        strNam2 = Base.get(d2,"Name","Error: no name") * strExt
    end
    
    f1 = FITSIO.FITS(filnam)
    m1 = FITSIO.read_header(f1[1])
    d1 = Base.read(f1[1])  # read an image from disk
    
    f2 = FITSIO.FITS(filnam2,"w")
    d2 = Base.write(f2, d1; header=m1)
    
    Base.close(f1)
    Base.close(f2)
        
    return filnam * " was saved as " * filnam2
end


"""
    fits_info(filnam [; info=false])

Metainformation of 'filnam.fits'
#### Example:
```
fits_info("T01.fits"; info=false)
T01.fits: file was found (for more information set info=true)
```
"""
function fits_info(filnam::String; info=false)
    
    _file_exists(filnam) || return "$filnam: file not found"

    f = FITSIO.FITS(filnam)
    m = FITSIO.read_header(f[1])
    d = Base.read(f[1])  # read an image from disk
    if info
        println(f)
        println("\r\n", f[1])
        Base.close(f)
        println("\r\nmetaInformation:\r\n", m)
    else
        Base.close(f)
        return filnam * ": file was found (for more information set info=true)"
    end
end

"""
    fits_key_create(filnam, key, value, comment)

create FITS key record
#### Example:
```
filnam = "T01.fits"
fits_key_create(filnam,"Test",3,"this is a test")
("TEST", 3, "this is a test")
```
"""
function fits_key_create(filnam::String, key::String, value, comment::String)
    
    key = Base.Unicode.uppercase(strip(key))
    
    _file_exists(filnam) || return "$filnam: file not found"

    _key_size_ok(key) || return error("error: $keynew: key name exceeds 8 characters (FITS record standard)")
    
    _key_available(filnam,key) || return "key '$key' cannot be created (already exists)"
    
    _key_comment_ok(comment) || println("Warning: comment truncated at 47 characters (FITS record standard)")
    
    f = FITSIO.FITS(filnam,"r+")                                    # open file append mode)
    k = FITSIO.write_key(f[1], key, value, Base.first(comment,47))  # write new key
    m = FITSIO.read_header(f[1])                                    # read hdu header from file   
    k = FITSIO.read_key(f[1],m.map[key])                            # test read new key
    
    Base.close(f)

    return k[1], k[2], k[3]
    
end

"""
    fits_key_delete(filnam, key)

edit FITS key record
#### Example:
```
filnam = "T01.fits"
fits_key_info(filnam,"test")
"key 'TEST' not in use"

fits_key_create(filnam,"Test",3,"this is a test")
("TEST", 3, "this is a test")

fits_key_delete(filnam,"test")
"TEST: key succesfully deleted"

fits_key_info(filnam,"test1")
"key 'TEST' not in use"
```
"""    
function fits_key_delete(filnam::String, key::String)
    
    key = Base.Unicode.uppercase(strip(key))
    
    _file_exists(filnam) || return "$filnam: file not found"
    
    _key_exists(filnam,key) || return "key '$key' cannot be deleted (key not found)"
    
    _key_reserved(filnam,key) && return error("error: $key: reserved key (FITS record standard)")
    
    buffer = "kanweg.fits"
    Base.Filesystem.isfile(buffer) && Base.Filesystem.rm(buffer)
    
    f1 = FITSIO.FITS(filnam,"r+")             # open file (read-only mode)
    d1 = Base.read(f1[1])                    # read data first hdu from file
    m1 = FITSIO.read_header(f1[1])           # read hdu header from file
    r1 = FITSIO.reserved_key_indices(m1)   
    
    f2 = FITSIO.FITS(buffer,"w")             # open file (write mode)
    d2 = Base.write(f2, d1)                  # write data first hdu to buffer (create default header)
    
    for i ∈ eachindex(m1.keys)
        k = FITSIO.read_key(f1[1],i)
        if i ∈ r1
            continue
        elseif k[1] == key
            continue
        else
            FITSIO.write_key(f2[1], k[1], k[2], k[3])
        end
    end
    
    Base.close(f2)    
    Base.close(f1)
    
    fits_copy(buffer,filnam)
    Base.Filesystem.rm(buffer)
    
    return _key_exists(filnam,key) || "$key: key succesfully deleted"
    
end


"""
    fits_key_edit(filnam, key, value, comment)

edit FITS key record
#### Example:
```
filnam = "T01.fits"
fits_key_edit(filnam, "USERTXT1", 3, "my text")
("USERTXT1", 3, "my text")
```
"""
function fits_key_edit(filnam::String, key::String, value, comment::String)
    
    key = Base.Unicode.uppercase(strip(key))
    
    _file_exists(filnam) || return "$filnam: file not found"
    
    _key_exists(filnam,key) || return "key '$key' not found"
    
    _key_reserved(filnam,key) && return error("error: $key: reserved key (FITS record standard)")
    
    _key_comment_ok(comment) || println("Warning: comment truncated at 47 characters (FITS record standard)")
    
    f = FITSIO.FITS(filnam,"r+")                                    # open file read-write mode)
    k = FITSIO.write_key(f[1], key, value, Base.first(comment,47))  # write new key
    m = FITSIO.read_header(f[1])                                    # read hdu header from file   
    k = FITSIO.read_key(f[1],m.map[key])                            # read new key
    
    Base.close(f)

    return k[1], k[2], k[3]
    
end

"""
    fits_key_info(filnam::String, key::String)

show FITS key record
#### Example:
```
filnam = "T01.fits"
fits_key_info(filnam,"NAXIS1")
("NAXIS1", 512, "length of data axis 1")
```
"""
function fits_key_info(filnam::String, key::String)
    
    key = Base.Unicode.uppercase(strip(key))
    
    _file_exists(filnam) || return error("error: $filnam: file not found")
    
    _key_exists(filnam,key) || return "key '$key' not in use"
    
    f = FITSIO.FITS(filnam,"r")                                     # open file read-write mode)
    m = FITSIO.read_header(f[1])                                    # read hdu header from file   
    k = FITSIO.read_key(f[1],m.map[key])                            # read new key
    
    Base.close(f)

    return k[1], k[2], k[3]
    
end

"""
    fits_key_rename(filnam, key, keynew)

rename FITS key
#### Examples:
```
filnam = "T01.fits"
fits_key_info(filnam,"test1")
"key 'TEST1' not in use"

fits_key_create(filnam,"Test",3,"this is a test")
("TEST", 3, "this is a test")

fits_key_rename(filnam,"test","test1")
("TEST1", 3, "this is a test")

fits_key_info(filnam,"test1")
("TEST1", 3, "this is a test")
```
"""    
function fits_key_rename(filnam::String, key::String, keynew::String)
    
    key = Base.Unicode.uppercase(strip(key))
    
    keynew = Base.Unicode.uppercase(strip(keynew))
    
    _file_exists(filnam) || return "$filnam: file not found"

    _key_size_ok(keynew) || return error("error: $keynew: key name exceeds 8 characters")
    
    _key_exists(filnam,key) || return "key '$key' cannot be renamed (key not found)"
    
    _key_reserved(filnam,key) && return error("error: $key: reserved key (FITS record standard)")
    
    (_key_available(filnam,keynew) || keynew == key) || return "$keynew: key already in use"
    
    buffer = "kanweg.fits"
    Base.Filesystem.isfile(buffer) && Base.Filesystem.rm(buffer)
    
    f1 = FITSIO.FITS(filnam,"r+")             # open file (read-write mode)
    d1 = Base.read(f1[1])                    # read data first hdu from file
    m1 = FITSIO.read_header(f1[1])           # read hdu header from file
    r1 = FITSIO.reserved_key_indices(m1)   
    
    f2 = FITSIO.FITS(buffer,"w")             # open file (write mode)
    d2 = Base.write(f2, d1)                  # write data first hdu to buffer (create default header)
    
    for i ∈ eachindex(m1.keys)
        k = FITSIO.read_key(f1[1],i)
        if i ∈ r1
            continue
        elseif k[1] == key
            FITSIO.write_key(f2[1], keynew, k[2], k[3]) 
        else
            FITSIO.write_key(f2[1], k[1], k[2], k[3])
        end
    end
    
    Base.close(f2)    
    Base.close(f1)
    
    fits_copy(buffer,filnam; protect=false)
    Base.Filesystem.rm(buffer)
    
    f = FITSIO.FITS(filnam,"r+")             # open file (read-write mode)
    m = FITSIO.read_header(f[1])             # read hdu header from file    
    k = FITSIO.read_key(f[1],m.map[keynew])  # read key
    
    Base.close(f)

    return k[1], k[2], k[3]
    
end
