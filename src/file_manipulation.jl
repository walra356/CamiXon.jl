using FITSIO

# A FITS file consists of a name and one or more header-data-units (HDUs) concatenated one after the other. The FITS object is a collection of these HDUs.

# The file 'filnamA.fits' is opened for 'reading only' by the command fileA = FITS(filnamA,'r').
# The second argument can be 'r' (read-only; default), 'r+' (read-write) or 'w' (read/create-write). 
# In the 'write' mode, a  new file is created or any existing file of the same name is overwritten. 

# The file fileA consists of one or more HDU's, fileA[1](, fileA[2],...). 
# HDUs consist of data blocks with a header containing metainformation, accessible by the commands read() and read_header().
# In the given example this becomes: dataA = read(fileA[1]) and metaInfoA = read_header(fileA[1])"
# The 3D array dataA represents a stack of images. The first image is obtained by the command imgA = dataA[:, :,1].

"""
    decompose_filnam(str)

Decompose filename into its name (and, if present, extension, prefix and numerator).
#### Examples:
```
strExample = "T23.01.fits"

dict = decompose_filnam(strExample)
Dict{String,String} with 4 entries:
  "Extension" => ".FITS"
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

    str = uppercase(str)
    ne = find_last(str,'.')[1]                 # ne: first digit of extension
    nl = length(str)                           # ne: length of file name including extension
    ne == nothing ? extension = false : extension = true

    if extension
        strNam = str[1:ne-1]
        strExt = str[ne:nl]
        isnumeric(str[ne-1]) ? n = ne-1 : n = nothing
        o = [("Name", strNam), ("Extension", strExt)]
    else
        strNam = str[1:ne-1]
        strExt = nothing
        isnumeric(str[nl]) ? n = nl : n = nothing
        o = [("Name", strNam)]
    end

    if n != nothing
        strNum = ""
        while isnumeric(str[n])
            strNum = str[n] * strNum
            n -= 1
        end
        strPre = str[1:n]
        Base.append!(o,[("Prefix", strPre),("Numerator", strNum)])
    end

    return Dict(o)

end

"""
    fits_combine(str1, str2; info=false)

Combine a series of .fits files into a single .fits file.
#### Example:
```
fits_combine("T01.fits", "T22.fits"; info=false)
T01-T22.FITS: file was created (for more information set info=true)
```
"""
function fits_combine(filnamFirst::String, filnamLast::String; info=false)

    dir = uppercase.(readdir())
    filnamFirst = uppercase(filnamFirst)
    filnamLast = uppercase(filnamLast)

    if filnamFirst ∉ dir
        error(filnamFirst * ": file not found")
    else
        d = decompose_filnam(filnamFirst)
        strPre = get(d,"Prefix","Error: no prefix")
        strNum = get(d,"Numerator","Error: no Numerator")
        strExt = get(d,"Extension","Error: no extension")
        valNum = parse(Int,strNum )
        numLeadingZeros = length(strNum) - length(string(valNum))
    end

    if filnamLast ∉ dir
         error(filnamLast * ": file not found")
    else
        d = decompose_filnam(filnamLast)
        strPre2 = get(d,"Prefix","Error: no prefix")
        strNum2 = get(d,"Numerator","Error: no Numerator")
        strExt2 = get(d,"Extension","Error: no extension")
        valNum2 = parse(Int,strNum2 )
        numLeadingZeros2 = length(strNum2) - length(string(valNum2))
    end

    if strPre ≠ strPre2
        error(strPre * " ≠ " * strPre2 * " (prefixes must be identical)")
    elseif strExt ≠ strExt2
        error(strExt * " ≠ " * strExt2 * " (file extensions must be identical)")
    elseif strExt ≠ ".FITS"
        error("file extension must be '.fits'")
    end

    numFiles = 1 + valNum2 - valNum
    fileFirst = FITS(filnamFirst)
    metaInfo = read_header(fileFirst[1])
    dataFirst = read(fileFirst[1])  # read an image from disk
    close(fileFirst)
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
        dataNext = read(fileNext[1])  # read an image from disk
        close(fileNext)
        dataStack[:, :,i] = dataNext[:, :,1]
    end

    filnamOut = strPre * strNum * "-" * strPre * strNum2 * strExt
    fileOut = FITS(filnamOut,"w")
    write(fileOut, dataStack; header=metaInfo)
    if info
        println("Output fileOut:\r\n", fileOut)
        println("\r\nOutput fileOut[1]:\r\n", fileOut[1])
        close(fileOut)
        println("\r\nmetaInformation:\r\n", metaInfo)
    else
        close(fileOut)
        return filnamOut * ": file was created (for more information set info=true)"
    end
end

"""
    fits_info(filnam; info=false)

Metainformation of 'filnam.fits'
#### Example:
```
fits_info("T01.fits"; info=false)
T01.fits: file was found (for more information set info=true)
```
"""
function fits_info(filnam::String; info=false)
    
    Base.Filesystem.isfile(filnam) ? true : return filnam * ": file not found in current directory"

    file = FITS(filnam)
    meta = read_header(file[1])
    data = read(file[1])  # read an image from disk
    if info
        println(file)
        println("\r\n", file[1])
        close(file)
        println("\r\nmetaInformation:\r\n", metaInfo)
    else
        close(file)
        return filnam * ": file was found (for more information set info=true)"
    end
end

"""
    fits_copy(filnam, filnamOut="")

Copy "filnam.fits" to "filnamOut.fits"
#### Examples:
```
fits_copy("T01.fits")
T01.fits was saved as T01 - Copy.fits

fits_copy("T01.fits","T01a.fits")
T01.fits was saved as T01a.fits
```
"""
function fits_copy(filnam, filnam2="")
    
    Base.Filesystem.isfile(filnam) ? true : return filnam * ": file not found in current directory"

    d1 = decompose_filnam(filnam)
    
    strNam = Base.get(d1,"Name","Error: no name")
    strExt = Base.get(d1,"Extension","Error: no extension")
        
    if filnam2 == ""
        filnam2 = strNam * " - Copy" * strExt
    else
        d2 = decompose_filnam(filnam2)
        strNam2 = Base.get(d2,"Name","Error: no name") * strExt
    end
    
    file1 = FITS(filnam)
    meta1 = read_header(file1[1])
    data1 = FITSIO.read(file1[1])  # read an image from disk
    
    file2 = FITS(filnam2,"w")
    dummy = FITSIO.write(file2, data1; header=meta1)
    
    close(file1)
    close(file2)
        
    return filnam * " was saved as " * filnam2
end
