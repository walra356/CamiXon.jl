module CamiXon

export decompose_filnam
export combine_fits_files
export find_all
export find_first
export find_last
export canonical_partitions
export integer_partitions

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
    combine_fits_files(str1, str2; info=false)

Combine a series of .fits files into a single .fits file.

#### Example:
```
combine_fits_files("T01.fits", "T22.fits"; info=false)
T01-T22.FITS: file was created (for more information set info=true)
```
"""
function combine_fits_files(filnamFirst::String, filnamLast::String; info=false)
    
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
        return println("jwError: " * strPre * " ≠ " * strPre2 * " (prefixes must be identical)")
    elseif strExt ≠ strExt2
        return println("jwError: " * strExt * " ≠ " * strExt2 * " (file extensions must be identical)")
    elseif strExt ≠ ".FITS"
        return println("jwError: file extension must be '.fits'")
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
        println(filnamOut * ": file was created (for more information set info=true)")
    end
end

"""
    find_all(A [,a...]; count=false)

A: string/array of elements of the same type

default   : Array containing the index (indices) of selected elements of A (default: all elements) 

count=true: The number of indices found for selected elements of A (default: all elements)  

#### Examples:
```
A = [:📑,:📌,:📢,:📌,:📞]
B = [1,2,3,2,5]
str = "aβcβd";

find_all(A) == find_all(B) == find_all(str)
true

find_all(A,:📌)
1-element Array{Array{Int64,1},1}:
 [2, 4]

find_all(str)
4-element Array{Array{Int64,1},1}:
 [1]
 [2, 4]
 [3]
 [5]

find_all(A; count=true)
4-element Array{Int64,1}:
 1
 2
 1
 1

str = "📑📌📢📌📞"
find_all(str,'📌')
1-element Array{Array{Int64,1},1}:
 [2, 4]
```
"""
function find_all(A::Union{String,AbstractArray{T,1}}, a::T...; count=false)  where T
    
    typeof(A) == String ? A = Base.collect(A) : nothing
    
    a == () ? a = Base.unique(A) : nothing
    
    o = [Base.findall(A .== Base.fill(a[i],length(A))) for i in eachindex(a)]
    
    return count ? length.(o) : o
    
end


"""
    find_first(A [,a...]; dict=false)

The first index of selected Array element

A: string/array of elements of the same type

default  : Array containing the first index (indices) of selected elements of A (default: all elements) 

dict=true: Dict for the first index (indices) of selected elements of A (default: all elements) 

#### Examples:
```
A = [:📑,:📌,:📢,:📌,:📞]
B = [1,2,3,2,5]
str = "aβcβd";

find_first(A) == find_first(B) == find_first(str)
true

find_first(A,:📌)
1-element Array{Array{Int64,1},1}:
 2

find_last(A,:📌; dict=true)
1-element Array{Pair{Symbol,Int64},1}:
 :📌 => 2

find_last(A; dict=true)
4-element Array{Pair{Symbol,Int64},1}:
 :📑 => 1
 :📌 => 2
 :📢 => 3
 :📞 => 5

find_first(str)
4-element Array{Int64,1}:
 1
 2
 3
 5
```
"""
function find_first(A::Union{String,AbstractArray{T,1}}, a::T...; dict=false)  where T
    
    typeof(A) == String ? A = Base.collect(A) : nothing
    
    a == () ? a = Base.unique(A) : nothing
    
    o = [Base.findfirst(A .== Base.fill(a[i],length(A))) for i in eachindex(a)]
    
    return dict ? [a[i] => o[i] for i in eachindex(a)] : o
    
end

"""
    find_last(A [,a...]; dict=false)

The last index of selected Array element

A: string/array of elements of the same type

default  : Array containing the lasst index (indices) of selected elements of A (default: all elements) 

dict=true: Dict for the lasst index (indices) of selected elements of A (default: all elements)

#### Examples:
```
A = [:📑,:📌,:📢,:📌,:📞]
B = [1,2,3,2,5]
str = "aβcβd";

find_last(A) == find_first(B) == find_first(str)
true

find_last(A,:📌)
1-element Array{Array{Int64,1},1}:
 4

find_last(A,:📌; dict=true)
1-element Array{Pair{Symbol,Int64},1}:
 :📌 => 4

find_last(A; dict=true)
4-element Array{Pair{Symbol,Int64},1}:
 :📑 => 1
 :📌 => 4
 :📢 => 3
 :📞 => 5

find_last(str)
4-element Array{Int64,1}:
 1
 4
 3
 5
```
"""
function find_last(A::Union{String,AbstractArray{T,1}}, a::T...; dict=false)  where T
    
    typeof(A) == String ? A = Base.collect(A) : nothing
    
    a == () ? a = Base.unique(A) : nothing
    
    o = [Base.findlast(A .== Base.fill(a[i],length(A))) for i ∈ eachindex(a)]
    
    return dict ? [a[i] => o[i] for i ∈ eachindex(a)] : o
    
end


function _canonical_partition(n::Int, m::Int)
    
    o = Base.fill(m,Base.cld(n,m))                              # init partition
    o[Base.cld(n,m)]=((n%m)≠0 ? n%m : m)                        # adjust last element of partition 
    
    return o 

end

"""
    canonical_partitions(n; header=false, reverse=true)

The canonical partition in integers of the integer n

header=true : unit patition included in output

#### Examples:
```
canonical_partitions(6; header=true, reverse=false)
6-element Array{Array{Int64,1},1}:
 [6]
 [5, 1]
 [4, 2]
 [3, 3]
 [2, 2, 2]
 [1, 1, 1, 1, 1, 1]

canonical_partitions(6; header=true)
6-element Array{Array{Int64,1},1}:
 [1, 1, 1, 1, 1, 1]
 [2, 2, 2]
 [3, 3]
 [4, 2]
 [5, 1]
 [6]

canonical_partitions(6)
5-element Array{Array{Int64,1},1}:
 [1, 1, 1, 1, 1, 1]
 [2, 2, 2]
 [3, 3]
 [4, 2]
 [5, 1]
```
"""
function canonical_partitions(n::Int, m=0; header=true, reverse=true)
    
    h = header ? n : n-1
    
    if m == 0
        if reverse
            o = [_canonical_partition(n,m) for m=1:h]
        else
            o = [_canonical_partition(n,m) for m=h:-1:1]
        end
    elseif 0 < m <= n
        o = _canonical_partition(n,m)
    else
        o = nothing
    end
    
    return o 

end



function _partition_count(n::Int,k::Int) 
    
    (n<0)|(k<0)|(k>n) ? 0 : (k==n)|(k==1) ? 1 : _partition_count(n-k,k) + _partition_count(n-1,k-1)
    
end

function _partition(a::Array{Int,1}, n::Int, i::Int, cp::Array{Array{Array{Int,1},1},1})
    
    o = a[1:i-1]
    m = a[i]-1                                           # m: partition value
    ni = n - Base.sum(o)                                 # ni: sub-partition index at partition index i 
      
    Base.append!(o,cp[ni][m])                            # complete partition by appending it to a
    
    return o
    
end

function _restricted_partitions(o::Array{Int,1}, n::Int, np::Int, ll::Array{Array{Int,1},1}, cp::Array{Array{Array{Int,1},1},1})
    
    oo = [o]
               
    for p=1:np-1
        i = Base.findlast(oo[p].>ll[length(oo[p])])
        Base.append!(oo,[_partition(oo[p],n,i,cp)])
    end
    
    return oo
    
end

"""
    integer_partitions(n [,m]; transpose=false, count=false)

default              : The integer partitions of n 

count=true           : The number of integer partitions of n 

transpose=false/true : for m>0 restricted to partitions with maximum part/length m

definitions:

The integer partition of the positive integer n is a nonincreasing sequence of positive integers p1, p2,... pk whose sum is n.

The elements of the sequence are called the parts of the partition. 

#### Examples:
```
integer_partitions(7)
15-element Array{Array{Int64,1},1}:
 [1, 1, 1, 1, 1, 1, 1]
 [2, 2, 2, 1]
 [3, 3, 1]
 [4, 3]
 [5, 2]
 [6, 1]
 [7]
 [2, 2, 1, 1, 1]
 [3, 2, 2]
 [4, 2, 1]
 [5, 1, 1]
 [2, 1, 1, 1, 1, 1]
 [3, 2, 1, 1]
 [4, 1, 1, 1]
 [3, 1, 1, 1, 1]

integer_partitions(7; count=true)
15

integer_partitions(7,4; count=true)
3

integer_partitions(7,4)
3-element Array{Array{Int64,1},1}:
 [4, 3]
 [4, 2, 1]
 [4, 1, 1, 1]

integer_partitions(7,4; transpose=true)
3-element Array{Array{Int64,1},1}:
 [2, 2, 2, 1]
 [3, 2, 1, 1]
 [4, 1, 1, 1]

```
"""
function integer_partitions(n::Int, m=0; transpose=false, count=false)
    
    ll = [ones(Int,l) for l=1:n]     
    cp = [canonical_partitions(m) for m=1:n] 
    oo = [ll[n]]
    pc = [_partition_count(n,m)  for m=1:n]
    
    np = m > 0 ? pc[m] : sum(pc)
    
    if !count
        
        if m == 0
            o = [_restricted_partitions(cp[n][p],n,pc[p],ll,cp) for p=2:n]
            for p=1:n-1 append!(oo,o[p]) end
        else
            oo = _restricted_partitions(cp[n][m],n,pc[m],ll,cp)
        end
        
        if transpose               
            for p=1:np
                l = length(oo[p])
                s=max(oo[p][1],l)
                mat = zeros(Int,s,s)
                for j=1:l for i=1:oo[p][j] mat[i,j]=1 end end
                oo[p] = [sum(mat[i,:]) for i=1:oo[p][1]]
            end
        
        end 

    end
    
    return count ? np : oo
    
end

end

