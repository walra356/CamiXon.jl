using FITSIO

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
        append!(o,[("Prefix", strPre),("Numerator", strNum)])
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
        println("jwError: " * filnamFirst * " (file not found)")
    else
        d = decompose_filnam(filnamFirst)
        strPre = get(d,"Prefix","Error: no prefix")
        strNum = get(d,"Numerator","Error: no Numerator")
        strExt = get(d,"Extension","Error: no extension")
        valNum = parse(Int,strNum )
        numLeadingZeros = length(strNum) - length(string(valNum))
    end

    if filnamLast ∉ dir
        println("jwError: " * filnamLast * " (file not found)")
    else
        d = decompose_filnam(filnamLast)
        strPre2 = get(d,"Prefix","Error: no prefix")
        strNum2 = get(d,"Numerator","Error: no Numerator")
        strExt2 = get(d,"Extension","Error: no extension")
        valNum2 = parse(Int,strNum2 )
        numLeadingZeros2 = length(strNum2) - length(string(valNum2))
    end

    if strPre ≠ strPre2
        return "Error: " * strPre * " ≠ " * strPre2 * " (prefixes must be identical)"
    elseif strExt ≠ strExt2
        return "Error: " * strExt * " ≠ " * strExt2 * " (file extensions must be identical)"
    elseif strExt ≠ ".FITS"
        return "Error: file extension must be '.fits'"
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
        L = length(filnamNext)
        filnamNext = strPre * "0"^numLeadingZeros * string(i) * ".fits"
        if L < length(filnamNext)
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

The metainformation of the file "filnam.fits"
#### Example:
```
fits_info("T01.fits"; info=false)         # assumption: file "T01.fits" in current directory
T01.FITS: file was found (for more information set info=true)
```
"""
function fits_info(filnam::String; info=false)
    
    dir = uppercase.(readdir())
    filnam = uppercase(filnam)
    
    if filnam ∉ dir
        return "Error: " * filnam * " (file not found)"
    else
        file = FITS(filnam)
        metaInfo = read_header(file[1])
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
end

"""
    fits_copy(filnam, filnamOut="")

Copy the file "filnam.fits"
#### Example:
```
fits_copy("T01.fits")                   # assumption: file "T01.fits" in current directory
T01.FITS was saved as T01 - Copy.FITS

fits_copy("T01.fits","T01a.fits")       # assumption: file "T01.fits" in current directory
T01.FITS was saved as T01A.FITS
```
"""
function fits_copy(filnam, filnamOut="")
    
    dir = uppercase.(readdir())
    filnam = uppercase(filnam)
    
    if filnam ∉ dir
        return println("Error: " * filnam * " (file not found)")
    else
        d = decompose_filnam(filnam)
        strNam = get(d,"Name","Error: no name")
        strExt = get(d,"Extension","Error: no extension")
        if filnamOut == ""
            filnamOut = strNam * " - Copy" * strExt
        else
            d2 = decompose_filnam(filnamOut)
            strNamOut = get(d2,"Name","Error: no name")
            filnamOut = strNamOut * strExt
        end
        file = FITS(filnam)
        metaInfo = read_header(file[1])
        data = read(file[1])  # read an image from disk
        fileOut = FITS(filnamOut,"w")
        write(fileOut, data; header=metaInfo)
        close(file)
        close(fileOut)
        return filnam * " was saved as " * filnamOut
    end
end
